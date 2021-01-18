import socket
import sys
import struct
import time
import datetime

##########
#


# dump: for debugging

def dump(s) :
    print ("DUMP:")
    for idx, x in enumerate(s) :
        if x >= ' ' and x <= '~' :
            c = x
        else :
            c = ' '
        print ("%2d %02x %c"%(idx,ord(x),c))

# block_parser: an object for parsing a response from the server.
# initialize it with a string (containing the server response) and
# the database connection object.  Then call getc(), get_n(), etc

class block_parser(object):

    def __init__(self, s, db ) :
        self.s = s
        self.db = db
        self.byte_order = db.byte_order
        self.cursor = 0

    # get a single byte as a string
    def getc( self ) :
        n = self.cursor
        if n >= len(self.s) :
            raise EOFError
        self.cursor += 1
        s = self.s[n]
        # print ("GETC",)
        # dump(s)
        # print ("cursor now",self.cursor)
        return s

    # get a fixed number of bytes as a string
    def get_n(self, count) :
        n = self.cursor
        if n+count >= len(self.s) :
            raise EOFError
        self.cursor += count
        s = self.s[n:n+count]
        # print ("GET_N",)
        # dump(s)
        # print ("cursor now",self.cursor)
        return s

    # get a null-terminated string
    def getstring( self ) :
        n = self.s[self.cursor:].find('\0')
        r = self.s[self.cursor:self.cursor+n]
        self.cursor += n
        # print ("GETSTRING",r)
        # print ("cursor now",self.cursor)
        return r

    # get an object header - returns tuple of ( classe, name )
    # classe is the integer class number 
    # name is a string
    def parse_object_header( self ) :
        while 1 :
            c = self.getc()
            if c == 'N' :
                classe = self.getc()
                name = self.getstring()
                return ( ord(classe), name )
            elif c == '/' :
                while c != '\n' :
                    c = self.getc()
            elif c == '\n' :
                pass
            elif c == '\0' :
                pass
            else :
                raise IOError

    def parse_table( self, table ) :
        # show -C begins it's output assuming that your cursor is at
        # the correct row, but column == -1.  The first time you
        # call parse_table, there are no rows or columns, so you start
        # at (0, -1).  On successive calls (like you might have with
        # tablemaker encores), you begin at the max row of the table, 
        # but still at column -1.

        row = len(table.rows)
        col = -1
        in_hash = 0

        while 1 :
            try :
                # print ("KTYPE",)
                ktype = self.getc()
            except EOFError :
                return

            if ktype == '#' :
                # end of record
                return

            elif ktype == '!' :
                #
                in_hash += 1

            elif ktype == '@' :
                #
                in_hash -= 1

            elif ktype == 'g' :
                # tag - 4 bytes
                value = self.get_n(4)
                value = struct.unpack(self.db.byte_order_u, value) [0]
                # print ("TAG", value, self.db.tag_num_to_name[value])
                value = self.db.tag_num_to_name[value]
                table.set( row, col, 'tag', value )

            elif ktype == 'K' or ktype == 'k' :
                # object - one byte class, string object name
                classe = ord(self.getc())
                name = self.getstring()
                # print ("OBJECT",classe, name)
                # BUG: not a nice representation
                table.set( row, col, 'object', ( classe, name ) )

            elif ktype == 'i' :
                # signed integer - 4 bytes
                value = self.get_n(4)
                value = struct.unpack(self.db.byte_order_i, value) [0]
                table.set( row, col, 'int', value )

            elif ktype == 'f' :
                # 4 byte float
                value = self.get_n(4)
                f = struct.unpack(self.db.byte_order_f, value) [0]
                table.set( row, col, 'float', f )

            elif ktype == 'd' :
                # 4 byte date - time_t
                value = self.get_n(4)
                f = struct.unpack(self.db.byte_order_u, value) [0]
                f = datetime.datetime.fromtimestamp( f )
                table.set( row, col, 'date', f )

            elif ktype == 'v' :
                # void : vNULL
                raise 'arf'     # not written yet

            elif ktype == 't' :
                # string
                value = self.getstring()
                # print ("STRING", value)
                table.set( row, col, 'string', value )

            elif ktype == '>' :
                # move to next column
                col += 1

            elif ktype == '.' :
                # move to next row
                table.copydown( row, row+1 )
                row += 1

            elif ktype == '<' :
                # increment row, decrement col
                table.copydown( row, row + 1 )
                row += 1
                col -= 1

            elif ktype == 'l' :
                # left: subtract %c from column.  Do not change row
                # this is always followed by a "." that increments the row 
                n = self.getc()
                col -= ord(n)

            elif ktype == '\n' :
                # newline is a nop
                pass

            elif ktype == '/' :
                # comment
                raise 'arf'     # not written yet

            elif ktype == 'D' or ktype == 'p' :
                # DNA
                # peptide
                # dna and peptide have the same format - a_data_len bytes of
                # unformatted data.
                if a_data_len == -1 :
                    raise Exception("server broken: sent type A data without length")
            
                raise 'arf'     # not written yet

            elif ktype == '\0' :
                # buggy servers sometimes insert a 00 byte
                pass

            else :
                # should not get here
                raise Exception("wac parse table not know ktype %d %c"%(ord(ktype),ktype))


##########
#

# this is an object that was found in the database.
# It has classe/name and a table of the data.  filled=True
# if we have the table, or False if we still need to
# get the table from the database.

class acedb_object(object) :

    def __init__( self, classe, name ) :
        self.classe = classe
        self.name = name
        self.filled = False
        self.table = None

    def show( self ) :
        print ("OBJECT",self.classe, self.name)
        self.table.show()

    def ac_has_tag( self, tagname ) :
        return self.table.find_tag( tagname )

    def ac_tag_int( self, tagname, default=None ) :
        return self.ac_tag_cell( tagname, 'int', default )

    def ac_tag_float( self, tagname, default=None ) :
        return self.ac_tag_cell( tagname, 'float', default )

    def ac_tag_text( self, tagname, default=None ) :
        return self.ac_tag_cell( tagname, 'string', default )

    def ac_tag_date( self, tagname, default=None ) :
        return self.ac_tag_cell( tagname, 'date', default )

    def ac_tag_cell( self, tagname, ctype=None, default=None ) :
        # find the item after the tag
        x = self.table.find_tag( tagname ) 
        if x is None :
            return default
        row, col = x
        c = self.table.cell( row, col+1 )
        if ctype is None :
            return c
        else :
            if c.dtype == ctype :
                return c.value
            else :
                return default

# This is the table that holds the object data.  Each cell of the
# table contains a acedb_table_cell

class acedb_table(object) :

    #
    def __init__( self, db=None ) :
        # current size
        self.cols = 0

        # array of arrays of cells
        self.rows = [ ]

        # the database we came from
        self.db = db

    #
    def show( self ):
        for n, row in enumerate(self.rows) :
            print ("Row",n)
            for m,v  in enumerate(row) :
                print ("    ",m,v.__str__())

    # ensure that there is space allocated for a particular row/column 
    def _grow_to( self, row, col ) :
        while row >= len(self.rows) :
            self.rows.append( [ ] )

        row_list = self.rows[row]
        while col >= len(row_list) :
            row_list.append( None )

    # get data from a specific row/col
    def cell( self, row, col, default=None ) :
        if row > len(self.rows) :
            return default
        if col > len(self.rows[row]) :
            return default
        return self.rows[row][col]

    # put data into a row/col
    def set( self, row, col, dtype, value ) :
        self._grow_to( row, col )
        self.rows[row][col] = acedb_table_cell( dtype, value )

    #
    def __str__( self ) :
        all_lines=[ ]
        for n in range(0,len(self.rows)) :
            all_lines.append( str(n)+ ":" + str(self.rows[n]) )
        return '\n'.join(all_lines)

    #
    def copydown( self, row_from, row_to ) :
        self._grow_to(row_to, 0)
        self.rows[row_to] = [ x for x in self.rows[row_from] ]

    # locate a tag in the table
    def find_tag( self, tagname ) :
        for rnum, r in enumerate(self.rows) :
            for cnum,c in enumerate(r) :
                # print (rnum, cnum, c.dtype, c.value)
                if c.dtype == 'tag' and c.value == tagname :
                    return rnum, cnum
        return None

# contents of a cell
#
# dtype is one of 'int', 'float', 'string', 'date', 

class acedb_table_cell(object) :
    def __init__( self, dtype, value ) :
        self.dtype = dtype
        self.value = value

    def __str__( self ) :
        return "ac_table_cell: %s - %s"%(self.dtype, self.value)

##########
#

# database access object
# db = acedb( 'a:localhost:12345::' )

class acedb(object) :

    def __init__( self, dbstr ) :
        '''dbstr is protocol:host:port:user:password
        '''
        l = dbstr.split(':')
        protocol = l[0]
        host=l[1]
        port=int(l[2])
        extras=l[3]
        if protocol == 'a' or protocol == 'acetcp' :
            self.db = acedb_connection_acetcp( protocol, host, port, extras)
        else :
            raise Exception('do not know specified protocol')

        self._ac_sys()

        self.db_date = 0

        # the object cache is indexed by the tuple (classe, name)
        self.object_cache = { }

    def ac_get_obj( self, classe, name, fillhint = False ) :
        # maybe we already have it
        idx = ( classe, name )
        if idx in self.object_cache :
            ob = self.object_cache[idx]
            if fillhint and not ob.filled :
                self.ac_fillobj( ob )
            return ob

        # look for it
        n = self._find_one( classe, name )
        if n == 0 :
            return None
        
        ob = acedb_object( classe, name )
        self.object_cache[idx] = ob
        if fillhint :
            self.ac_fillobj_part2( ob )

        return ob

    def ac_fillobj( self, ob ) :
        # fill the object, 
        n = self._find_one( classe, name )
        if n == 0 :
            return

        self._fillobj_part2( ob )

    def ac_fillobj_part2( self, ob ) :
        # ac_fillobj_part2 () assumes that the object is currently the only object
        # in the active keyset.  That is OK because either ac_fillobj () or ac_get ()
        # caused that to be true before calling here.
        #
        # This entry point only exists so we can eliminate the server transaction
        # in ac_fillobj () when you ac_get () an object.
        ob.table = acedb_table()

        s = self.ac_command( 'show -C' )
        parse = block_parser(s, self)

        # print ("PARSE")
        # dump(s)

        # must consume object header even though it contains nothing of interest
        classe, object_name = parse.parse_object_header()

        # AceC checks that the object name and class number are what we expect
        # AcePython does not

        ob.get_date = self.db_date
        self.db_date += 1

        # flag filled before we call into parse table
        ob.filled = True

        parse.parse_table( ob.table )


    def _find_one( self, classe, name ) :
        s = self.ac_command( 'query find %s IS %s'%(classe, name))
        return self.active_objects_count(s)

    def active_objects_count( self, s ) :
        try :
            s = s.split('Active Objects')
            s = s[0]
            s = s.split('//')
            s = s[-1]
            return int(s)
        except :
            raise Exception('server sent invalid response for active objects')

    def _ac_sys( self ) :
        b = self.ac_command("system_list")

        # pick out the byte order
        n = struct.unpack_from('<i',b)[0]
        if n == 0x12345678 :
            self.byte_order = '<'
        else :
            n = struct.unpack_from('>i',b)[0]
            if n == 0x12345678 :
                self.byte_order = '>'
            else :
                raise Exception('unrecognized byte order %x'%n)

        self.byte_order_f = self.byte_order + 'f'
        self.byte_order_u = self.byte_order + 'I'
        self.byte_order_i = self.byte_order + 'i'

        # skip over 
        #   4 - the byte order marker 
        #   3 - the number of tags
        #   1 - '>' at the beginning of the tags
        b = b[8:]

        # the result is line oriented
        b = b.split('\n')

        self.tag_num_to_name = { }
        self.tag_name_to_num = { }
        self.classe_num_to_name = { }
        self.classe_name_to_num = { }

        # tag #0 is not used; it also comes out in a different format,
        # so we don't recognize it.  So, start ntags=1
        ntags = 1

        # 
        nclasses = 0
        for line in b :
            if line.startswith('.g') :
                # load a tag - '.gXXX'+'\0'
                self.tag_num_to_name[ntags] = line[2:-1]
                self.tag_name_to_num[line[2:-1]] = ntags
                ntags += 1
            elif line.startswith('.k') :
                # load a class - '.kXXX'+'\0'
                # blank class name means a non-meaningful class number, so we don't
                # need to remember it.
                c = line[2:-1]
                if c != '' :
                    self.classe_num_to_name[nclasses] = line[2:-1]
                    self.classe_name_to_num[line[2:-1]] = nclasses
                nclasses += 1
            elif line.startswith('#') :
                # end of meaningful part of response
                break
            else :
                # the first line goes unrecognized
                pass

    def ac_command( self, cmd ) :
        msg = self.db.partial_command( cmd )
        l = [ msg ]
        while self.db.encore :
            msg = self.db.partial_command('encore')
            l.append( msg )
        return ''.join(l)


# generic class of a connection object - can be used with isinstance
class acedb_connection(object):
    pass

# the acetcp connection object - I don't know if we will ever
# implement the RPC or the Sanger server protocols because they are
# harder to do in python.  (And this works.)

class acedb_connection_acetcp(acedb_connection) :

    def __init__( self, protocol, host, port, extras, timeout = None ) :
        # pick out timeout even though we don't use it; maybe someday
        if '/' in protocol :
            self.timeout = int( protocol.split('/')[1] )
        else :
            self.timeout = 300

        # the database specifier is "acetcp:host:port" or "a:host:port"
        # but thanks to HTTP, people want to write "acetcp://host:port".
        # I will let them.
        if host.startswith('//') :
            host = host[2:]

        if host == '' :
            host = 'localhost'

        #
        self.socket = socket.create_connection( (host,port), timeout )

        self.socket.setsockopt( socket.SOL_SOCKET, socket.SO_SNDBUF, 49152 )
        self.socket.setsockopt( socket.SOL_SOCKET, socket.SO_RCVBUF, 49152 )

        # 6 is TCP
        self.socket.setsockopt( 6, socket.TCP_NODELAY, 1 )

        self.incoming_transaction = incoming_transaction(self.socket)
        self.outgoing_transaction = outgoing_transaction(self.socket)

        # consume the challenge that the server will send us on connect.  if
        # you implement authentication, store the challenge to use in later
        # packets.
        msg_type, msg = self.incoming_transaction.read( )
        
        if msg_type == 'C' :
            self.challenge = msg
        elif msg_type == 'E' :
            raise Exception(msg)
        else :
            raise Exception('server sent invalid response')

    def close( self ):
        'ac_close_acetcp'
        pass

    def partial_command( self, cmd ):

        # recognize the special command "encore" and make it ask for an
        # 'encore' transaction instead of a command transaction.  This is a
        # yucky way to recognize it, but it is compatible with all the other 
        # transports.
        if cmd == 'encore' :
            self.outgoing_transaction.send( 'e', cmd + '\0' )
        elif cmd == 'comment' :
            self.outgoing_transaction.send( 'n', cmd + '\0' )
        else :
            self.outgoing_transaction.send( 't', cmd + '\0' )

        # 
        self.encore = 0

        # read a response
        msg_type, msg = self.incoming_transaction.read()
        if msg_type == 'C' :
            raise Exception('server challenged us again')
        elif msg_type == 'a' :
            raise Exception('server selected authentication - not implemented')
        elif msg_type == 'E' :
            raise Exception('server sent back error - %s'%msg)
        elif msg_type == 'R' :
            self.encore = ord(msg[0])
            return msg[1:]
        else :
            raise Exception('server sent unrecognized transaction type %c %02x'%(msg_type,ord(msg_type)))

    def lazy_command( self, cmd ):
        # send a command the expects no response
        self.outgoing_transaction.send( 'n', cmd + '\0' )


# modeled after wtcp/tcp_transactions.c

class incoming_transaction(object) :
    BUFFER_SIZE = 49152

    def __init__( self, socket ) :
        self.eofcount = 0
        self.data = ''
        self.awaiting = 0
        self.socket = socket

    def read( self ) :
        t = self.read_op()
        while t is None :
            t = self.read_op()
            # bug: sleep?
        return t[0], t[1:]

    def read_op( self, callback=None, cookie=None, limit_one=True ) :
        '''
        Return value is:
            None - no data yet; try again
            str - a message
        on EOF, raises EOFError
        '''
        l = len(self.data)
        self.data = self.data + self.socket.recv( self.BUFFER_SIZE )
        if len(self.data) == l :
            self.eofcount += 1
            if self.eofcount > 3 :
                raise EOFError
            else :
                return None    # no data, not declared EOF yet

        while 1 :

           # a transaction has at least 3 bytes of length at the
           # begining
            if len(self.data) < 3 :
                return None

            # If self.awaiting is 0, we don't know how many bytes we are
            # expecting.  The length is at the front of the buffer, so
            # we can just pick it up.
            if self.awaiting == 0 :
                self.awaiting = ord(self.data[0]) + ( ord(self.data[1]) << 8 ) + ( ord(self.data[2]) << 16 )

            # if the size of the data we have is smaller than the length
            # of the transaction, we can't process it.  We return, and
            # when more data is available, we will be called again to 
            # again read from the buffer.
            if len(self.data) < self.awaiting :
                return None

            # Ok, there is a whole transaction ready.  Process it and then
            # tear it off the front of the buffer
            result = self.data[0:self.awaiting]
            self.data = self.data[self.awaiting:]
            self.awaiting = 0

            if callback is None :
                return result[3:]
            else :
                callback(result, cookie)
            if limit_one :
                return 0


class outgoing_transaction(object) :
    def __init__( self, socket ) :
        self.socket = socket
        self.auth = None

    def send( self, msg_type, message ) :
        # length is 3 bytes
        # msg_type is 1 byte
        total_len = len(message) + 3 + 1;

        # user can give message type as an int or single char string
        if not isinstance( msg_type, int) :
            msg_type = ord(msg_type)

        #
        msg_l = [ 
            chr(total_len & 0xff),
            chr((total_len >> 8) & 0xff),
            chr((total_len >> 16) & 0xff),
            chr(msg_type),
            message
            ]
        self.socket.send( ''.join(msg_l) )


if __name__ == '__main__' :
    db = acedb( 'a:localhost:12345::' )

    s = db.ac_get_obj( 'arf', 'self_ref', 1 )
    s.show()
    #
    s = db.ac_get_obj( 'arf', 'a', 1 )
    # print (s)

    # s = db.ac_command( 'list -C' )
    # print (s)



