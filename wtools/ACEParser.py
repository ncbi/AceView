#!/usr/bin/python

## Todo List: - complete the predefined classes (what about ?Bat etc. ?) (hopefully already done)
##            - store comments from model files to write them back out with the model
##            - add methods to use the model class to represent a java swing jtree
##            - add COORD modifier and @ ????
##            - The CLASSDICT should be created as the classes are created, maybe the classes should be stored in a dictionary instead of a list
##            - add the @ modifier
                
import string  ## imports a module for string operations
import sys     ## imports a module to get command line arguments


class Model:   ## This class represents the model as a tree    
    ## Systemmodels is a long string that hold the textual representation of the classes and datatypes defined in sysclasses.c
    ## This might be necessary for the validation.
    ## At least ?Colour and ?View are used in some models

    Systemmodels=["?Model","?Text","?Comment","?Keyword","?Class","?Tag","SourceCode","?Include","?Session","?Display","?Table","#Table_definition","?Colour","?UserSession","?Jade","?View"]

    
    CharacterArrays=["?LongText","?DNA","?Peptide","?BasePosition","?BaseQuality","?URL","?BaseCall","?KeySet","?MatchTable"]   ## some character arrays of ACEDB.
                      ## they are needed for the validation process

    BasicDataTypes=["Int","Text","Float","DateType"]
    ## those are needed to find out if an entry is a label or a datatype 
    
    def __init__(self,FileContent):  ## Constructor
        self.FileContent=FileContent ## store the filecontent so we can find linenumbers later on
        self.Classes=[]              ## A List that holds all the classes
        self.DefinedDataTypes=[]     ## A list that holds the user defined datatypes
        Blocks=self.MakeBlocks(FileContent)  ## Divides the content of the file into blocks which correspond to the individual classes
        for Block in Blocks: 
            if (Block[0][0]== "?"):   ## if it starts with a '?', it is a class
                self.Classes.append(ClassType(self,Block,self.BasicDataTypes)) ## add a new class instance to Classes
            else:
                if (Block[0][0]=="#"): ## if it starts with a '#' it is a datatype
                    self.DefinedDataTypes.append(UserDefinedType(self,Block,self.BasicDataTypes)) ## add a new datatype to the user defined datatypes
                else: ## if it doesn't start with '?' or '#' than it is an error
                    print "line "+str(self.FindLine(self.GetPath()))+"ERROR: "+self.GetPath()+" This Datatype is unknown: "+string.split(Block[0])[0]
                    print "Please check your model."
                    sys.exit(0) ## exit the program
        
    def MakeBlocks(self,ModelText):  
        ClassBlocks=[]
        ClassBlock=[]        
        Text=string.join(ModelText,"")
        while (string.find(Text,"/*") != -1): ## removes multi line comments 
            x=string.find(Text,"/*")
            y=string.find(Text,"*/")
            if (x > y): ## since it doesn't work for all cases print error message and exit
                print "There is some problem with the comments."
                print "Please remove all nestet comments."
                sys.exit(0)
            else:
                Text=Text[:x]+Text[y+2:]
        ModelText=string.split(Text,"\n")
        NumOfClasses=0
        for Line in ModelText:
            LINENUMBER=ModelText.index(Line)+1
            x=string.find(Line,"//")
            if (x != -1):
                Line=Line[:x]   ## removes one line comment
                Comment=Line[x:]
            if (len(string.split(Line)) > 0):
                Line=string.expandtabs(Line)
                if (Line[0] not in [" ","\t"]):  ## if the line starts with anything else than ' ' or '\t' a new class starts
                    if ((Line[0]) not in ["#","?"]):
                        print "line "+str(LINENUMBER)+": ERROR. This is neither a class nor a datatype."
                    else:
                        NumOfClasses=NumOfClasses+1
                        if (NumOfClasses !=1):    ## if it is not the first class
                            ClassBlocks.append(ClassBlock)  ## add the class
                        ClassBlock=[Line]  ## this Line belongs to the next class
                else:
                    ClassBlock.append(Line) ## it is the first time a `?' or `#' is encountered, so there is nothing to add yet
        ClassBlocks.append(ClassBlock)      ## We have to add the last class block before we return to the main program
        return ClassBlocks  ## this return a list of Textblocks that correspond to classes

    def __str__(self):  ## this represents the model if you say 'print 'instance of Model'
        String=""
        for Class in self.Classes:
            String=String+str(Class)
        for UDataType in self.DefinedDataTypes:
            String=String+str(UDataType)
        return String

    def __repr__(self):   ## similar to __str__
        return self.__str__()


    def toString(self):   ## actualy not used, but might come handy when interfacing to java
        return self

    def Validate(self):   ## this routine checks the model
        ERROR=""
        ValidIdentifiers=[]
        CLASSDICT=self.MakeClassDict()  ## this makes a dictionary of all the defined classnames and pointers to the corresponding class in the model
        for CLASS in self.Classes:  ## let every class validate itself
            ERROR=ERROR+CLASS.Validate(self.BasicDataTypes,CLASSDICT,self.CharacterArrays+self.Systemmodels)
        for DATATYPE in self.DefinedDataTypes: ## let every datatype validate itself
            ERROR=ERROR+DATATYPE.Validate(self.BasicDataTypes,CLASSDICT,self.CharacterArrays+self.Systemmodels)
        return ERROR
    
    def MakeClassDict(self):  ## this makes a dictionary of all classnames and pointers to the corresponding classes
                              ## with this it is possible to get a classinstance with the name of the class
                              ## {classname:pointer_to_class}
        CLASSDICT={}
        for CLASS in self.Classes+self.DefinedDataTypes:
            CLASSDICT[CLASS.GetName()]=CLASS
        return CLASSDICT

    def GetPath(self):  
        return ""

    def FindLine(self,ObjectPath):  ## find the linenumber of an object in the filecontent
        LINENUMBER=-2
        Object=string.split(ObjectPath," -> ")
        if (Object[0][0] not in ["#","?"]):
            LINENUMBER=-3
        else:
            START=-1
            END=-1
            for Line in self.FileContent:
                if ((Line[0] in ["#","?"]) and (string.split(Line)[0] == Object[0])):
                    START=self.FileContent.index(Line)
                    continue
                if ((Line[0] in ["#","?"]) and (START != -1)):
                    END=self.FileContent.index(Line)-1
                    break
            if (START != -1):
                if (END == -1):
                    END=len(self.FileContent)-1
                OBJECTCOUNT=1
                i=START
                while (i <= END):
                    if (string.find(string.expandtabs(self.FileContent[i])," "+Object[OBJECTCOUNT]) != -1): ## this is the most simple appoach, but errorprone
                        OBJECTCOUNT=OBJECTCOUNT+1    ## should be changed soon
                        i=i-1
                        if (OBJECTCOUNT==len(Object)):
                            LINENUMBER=i+2
                            break
                    i=i+1
                        
        return LINENUMBER
                    
                    
                    

class ClassType:  ## this class represents a model class
    def __init__(self,Parent,ClassBlock,BasicDataTypes): ## Constructor
        self.Parts=[]  ## Parts is a list that stores the branches (labels) and leafes (datatypes)
        self.SetParent(Parent) ## stores the parent node
        self.SetName(string.split(ClassBlock[0])[0])  ## sets the name of the class
        NextLevel=self.GetNextLevelBlocks(ClassBlock) ## makes the next level of hirearchy
             ## this just makes textblocks of all the entries that belong to the next level in the tree
        for Object in NextLevel: ## for each textblock try to make an object and add it to the list
            ## the entries are created as lists of lists
            ## each line in the modelfile corresponds to a list
            ## lists that correspond to datatypes will have more entries in one list
            ## list corresponding to labels have only one entry
            ## [
            ##  [label1],
            ##  [label2],
            ##  [datatype1, datatype2, datatype3...],
            ##  [label4]
            ## ]
            
            O=self.MakeObject(Object,BasicDataTypes)
            if (O):
                self.AddObject(O)
            else: ## if no object could be made, print an error and exit
                print "This block does not represent a valid entry:"
                print Object
                sys.exit(0)

    def GetNextLevelBlocks(self,Block):  ## makes textblocks corresponding to the next level in the tree
        BLOCKS=[] 
        TEMP=string.split(string.strip(Block[0]))  ## split first line into parts
        if (len(TEMP) > 1):                        ## if object has parts
            INDENT=string.find(Block[0]," "+TEMP[1])+1
            TEMP2=[]
            for LINE in Block:
                if (LINE[INDENT-1] not in [" ","\t"]):     ## if the character before the current indentation level is a non whitespace, then we have a problem with the indentation
                    print self.GetPath()+": There seems to be a problem with the indentation.\n\n"
                    sys.exit(0)
                else:
                    LINE=LINE[INDENT:]      ## shorten each line according to the indentation level
                    if (LINE[0] != " "):    ## if the line doesn't start with a space, we have the start of a new object  
                        if (TEMP2 != []):   ## is it the first one ?, if not
                            BLOCKS.append(TEMP2)    ## add the current object to the list
                        TEMP2=[LINE]
                    else:
                        TEMP2.append(LINE)     ## starts with a space? So it is part of the current object
            BLOCKS.append(TEMP2)               ## don't forget to add the last one to the list
        else:                                  ## if object doesn't have parts, it can only one line long
                                               ## otherwise there is a problem with the indentation
            if (len(Block)>1):
                print self.GetPath()+": There seems to be a problem with the indentation.\n\n"
                sys.exit(0)
        return BLOCKS

    def __repr__(self):  ## similar to __str__
        return self
    
    def __str__(self):   ## returns a string representation of the class
        STRING=""
        if (len(self.Parts) > 0): 
            for PART in self.Parts:  ## let every entry represent itself
                for P in PART:
                    STRING=STRING+str(P)+" "
                STRING=STRING[:-1]
            LINES=string.split(STRING,"\n")
            LINES=LINES[:-1]
            LINES[0]=self.GetName()+"\t"+LINES[0]  ## add the name of this instance to the beginning of the first line
            INDENT=(len(self.GetName())/8)+1     ## Indentation of the other lines
            for X in range(1,len(LINES)):
                LINES[X]=INDENT*"\t"+LINES[X]    ## indent all the other lines
            STRING=string.join(LINES,"\n")
            STRING=STRING+"\n"
        else:
            STRING=self.GetName()+"\n"           ## if the instance doesn't have entries, return just the name
        return STRING
   
    def SetParent(self,Parent):  ## sets the parent node
        self.Parent=Parent

    def GetParent(self):         ## return the parent node
        return self.Parent

    def SetParts(self,Parts):    ## not used, but sets Parts
        self.Parts=Parts

    def GetParts(self):
        return self.Parts        ## returns all entries

    def HasLabel(self,LabelName):  ## checks wether the class has a specified label
        LABEL=""
        for PART in self.Parts:
            if (string.upper(PART[0].GetName()) == string.upper(LabelName)):  ## if the class has the label
                LABEL=PART[0].GetPath()           ## return the Path to the label
                break
        if (LABEL == ""):                         ## if the class doesn't contain the label directly, check further down
            for PART in self.Parts:
                LABEL=PART[0].HasLabel(LabelName)
                if (LABEL != ""):
                    break
        return LABEL  ## returns an empty string if the label wasn't found

    def Validate(self,BasicDataTypes,ClassDict,CharacterArrays): ## checks the class
        ERROR=""
        for PART in self.Parts:  ## let every entry check itself
            for P in PART:
                ERROR=ERROR+P.Validate(BasicDataTypes,ClassDict,CharacterArrays)
        if (self.GetName()[0] != "?"): ## if the name doesn't start with an `?`, it's not a class
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". Classes have to start with an '?'.\n"        
        if (string.upper(self.GetName()) in BasicDataTypes+["UNIQUE","XREF","REPEAT"]):
            ## if the name corresponds to a reserved keyword, than print an error message
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". The name "+self.GetName()+" is not allowed for a class.\n"
        if (("#"+string.upper(self.GetName()[1:])) in map(string.upper,ClassDict.keys()+CharacterArrays)): ## checks wether there already exists a defined datatype with the same name
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". There already exists a defined datatype with this name.\n"
        if (len(self.Parts) < 1): ## if the class is empty
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". Classes must have labels. They can not be empty.\n"
        return ERROR

    def SetName(self,Name): 
        self.Name=Name

    def GetName(self):
        return self.Name

    def GetPath(self):      ## returns the path of the instance ?Classname -> Label1 -> Label11 -> Datatype111
        PATH=""
        if (hasattr(self,"Parent")):
            PATH=self.Parent.GetPath()
        if (PATH==""):
            PATH=self.Name
        else:
            PATH=PATH+" -> "+self.Name
        return PATH

    def AddObject(self,Object):    ## just adds an entry
        self.Parts.append(Object)

    def RemoveObject(self,Object): ## removes an entry
        self.Parts.remove(Object)

    def MakeDataPackages(self,Object,BasicDataTypes):  ## if the object is a line with several datatypes, it has to be divided into the individual datypes
        ## this makes lists of the datatypes and their modifiers
        BITS=string.split(Object[0])
        PACKAGES=[]
        LAST=-1
        for X in range(0,len(BITS)):
            if ((string.upper(BITS[X]) in map(string.upper,BasicDataTypes)) or (BITS[X][0] in ["#","?"])):
                if (LAST != -1):
                    if (string.upper(BITS[X-1]) == "UNIQUE"):
                        PACKAGE=BITS[LAST:X-1]
                        PACKAGES.append(PACKAGE)
                    else:
                        PACKAGE=BITS[LAST:X]
                        PACKAGES.append(PACKAGE)
                    LAST=X
                else:
                    LAST=0
        if (string.upper(BITS[LAST-1]) =="UNIQUE"):
            PACKAGES.append(BITS[LAST-1:])
        else:
            PACKAGES.append(BITS[LAST:])
        return PACKAGES
                                        
                                     

    def MakeObject(self,Object,BasicDataTypes):  ## this routine checks wether an object is a label or a datatype or a 'UNIQUE list'
        IDENTIFIER = string.split(Object[0])[0]
        OBJECT=[]
        if ((string.upper(IDENTIFIER) in map(string.upper,BasicDataTypes)) or (IDENTIFIER[0] in ["#","?"])):
            PACKAGES=self.MakeDataPackages(Object,BasicDataTypes)
            for PACKAGE in PACKAGES:
                OBJECT.append(DataType(self,PACKAGE,BasicDataTypes))
        else:
            if (string.upper(IDENTIFIER) == "UNIQUE"):
                if ((len(Object) == 1) and ((string.upper(string.split(Object[0])[1]) in map(string.upper,BasicDataTypes)) or (string.split(Object[0])[1][0] in ["#","?"]))):
                    PACKAGES=self.MakeDataPackages(Object,BasicDataTypes)
                    OBJECT=[]
                    for PACKAGE in PACKAGES:
                        OBJECT.append(DataType(self,PACKAGE,BasicDataTypes))
                else:
                    OBJECT.append(ExclusiveListType(self,Object,BasicDataTypes))
            else:
                if (string.upper(IDENTIFIER) in ["XREF","REPEAT"]):
                    print self.GetPath()
                    print "Fatal error. This entry can not have a modifier "+IDENTIFIER
                    sys.exit(0)
                else:
                    OBJECT.append(LabelType(self,Object,BasicDataTypes))
        return OBJECT

    def GetNumParts(self):    ## returns the number of entries
        return len(self.Parts)

    def GetLenPart(self,INDEX):            ## returns the length of an entry
                                           ## for labels this is 1
                                           ## for lines with datatypes, this corresponds to the number of datatypes on the line
        if (self.GetNumParts() >= INDEX):
            RESULT=len(self.Parts[INDEX])
        else:
            RESULT=-1
        return RESULT

    def FindElement(self,Element):     ## this searches for an instance of a label or datatype in Parts and returns its index
        RESULT=(-1,-1)
        for PART in self.Parts:
            if (Element in PART):
                RESULT=(self.Parts.index(PART),PART.index(Element))
        return RESULT


    def GetLineNumber(self,ClassDict):
        MODEL=ClassDict[string.split(self.GetPath())[0]].GetParent()    
        if (isinstance(MODEL,Model)):
            LINENUMBER=MODEL.FindLine(self.GetPath())
        else:
            LINENUMBER=-1
        return LINENUMBER
        


class LabelType(ClassType):  ## this class represents a label
                             ## it inherits everything from the ClassType class
                             ## only the validation method is overwritten
    def Validate(self,BasicDataTypes,ClassDict,CharacterArrays): ## validates the label
        ERROR=""
        for PART in self.Parts: ## let all entries check themselves
            for P in PART:  
                ERROR=ERROR+P.Validate(BasicDataTypes,ClassDict,CharacterArrays)
        if ((self.GetName()[0] in ["#","?"]) or (string.upper(self.GetName()) in map(string.upper,BasicDataTypes)+["UNIQUE","XREF","REPEAT"])):
            ## if the name starts with `#` or `?` or corresponds to a reserved keyword, than print an error message
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))++": ERROR  "+self.GetPath()+". This name is not allowed for a label.\n\n"
        return ERROR
    
class ExclusiveListType(ClassType):  ## this class represents a UNIQUE list of labels
                                     ## UNIQUE datatypes are just datatypes with the modifier UNIQUE
                                     ## this class also inherits everything from ClassType and overwrites the validation method
    def Validate(self,BasicDataTypes,ClassDict,CharacterArrays):
        ERROR=""
        for PART in self.Parts:
            if (len(PART) > 1): ## if the entry has more than one entry, it can not be a label
                ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". UNIQUE lists can only have labels not datatypes.\n"
            for P in PART:
                ERROR=ERROR+P.Validate(BasicDataTypes,ClassDict,CharacterArrays)
        if (string.upper(self.GetName()) != "UNIQUE"): ## actualy the UNIQUE list is just a label that has the name UNIQUE
                                         ## if it doesn't have that name, it's an error
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". UNIQUE lists must have name 'UNIQUE'.\n"
        if (len(self.Parts) < 1):  ## UNIQUE lists can not be empty
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". UNIQUE lists must have labels. They can not be empty.\n"
        return ERROR

class UserDefinedType(ClassType):  ## this class represents the user defined datatypes
                                   ## it inherits everything from ClassType and overwrites the validation method
    def Validate(self,BasicDataTypes,ClassDict,CharacterArrays):
        ERROR=""
        for PART in self.Parts: ## let every entry check itself
            for P in PART:
                ERROR=ERROR+P.Validate(BasicDataTypes,ClassDict,CharacterArrays)
        if (self.GetName()[0] != "#"): ## per definition user defined datatypes have to start with an '#'
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". User defined datatypes  have to start with an '#'.\n"
        if (len(self.Parts) < 1):  ## an empty user defined datatype doesn't make sense
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". User defined datatypes not be empty.\n"
        if ((string.upper(self.GetName()[1:]) in map(string.upper,BasicDataTypes)+["UNIQUE","XREF","REPEAT"])):
            ## if the name corrsponds to a reserved keyword, than print an error message
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR  "+self.GetPath()+". This name is not allowed for a user defined datatype.\n\n"
        if ("?"+string.upper(self.GetName()[1:]) in map(string.upper,ClassDict.keys()+CharacterArrays)): ## checks wether there already exists a class with the same name
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". There already exists a class with this name.\n"
        return ERROR
    


class DataType:  ## this class represents a datatype
    def __init__(self,Parent,Object,BasicDataTypes):  ## Constructor
        self.Modifiers=[]   ## this lists hold the modifiers
        self.SetParent(Parent)  ## sets the parent node
        for P in Object:        ## an initialisation a datatype object gets a list of string that represents the individual parts of the datatype, like modifiers, name etc.
                                ## this list is not ordered
            if (string.upper(P) == "UNIQUE"):
                self.AddModifier(UNIQUE_Modifier(self,P))
            if ((string.upper(P) in map(string.upper,BasicDataTypes)) or (P[0] in ["#","?"])):
                self.SetName(P)
            if (string.upper(P) == "XREF"):
                self.AddModifier(XREF_Modifier(self,P,Object[Object.index(P)+1]))
            if (string.upper(P) == "REPEAT"):
                self.AddModifier(REPEAT_Modifier(self,P))
            

    def SetParent(self,Parent):  ## sets the parent node
        self.Parent=Parent

    def GetParent(self):         ## returns the parent node 
        return self.Parent

    def Validate(self,BasicDataTypes,ClassDict,CharacterArrays):  ## checks the datatype
        ERROR=""
        if ((string.upper(self.GetName()) not in map(string.upper,BasicDataTypes)) and (self.GetName()[0] not in ["?","#"])):
            ## the name has to start with `?' or `#` or be a basic data type
            ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". This datatype is not valid.\n"
        else:
            if (self.GetName()[0] in ["#","?"]):  ## if it is a class, it must be defined in the model
                if ((string.upper(self.GetName()) not in map(string.upper,ClassDict.keys()+CharacterArrays)) and ("?"+string.upper(self.GetName()[1:]) not in map(string.upper,ClassDict.keys()+CharacterArrays))):
                    ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+". This class or datatype was never defined.\n"
                else: ## check if the classname has consistens spelling concerning upper and lower case letters
                    if ((self.GetName() not in ClassDict.keys()+CharacterArrays) and ("?"+self.GetName()[1:] not in ClassDict.keys()+CharacterArrays)):
                        ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": WARNING! "+self.GetPath()+": Please check the spelling.\n"
##                if (not isinstance(self.GetParent(),LabelType)):
##                  ERROR=ERROR+"ERROR in "+self.GetPath()+": Datatypes must follow a label not a classname.\n"
            else:  ## check if the basic datatype has consistent usage of upper and lower case letters
                if (self.GetName() not in BasicDataTypes):
                        ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": WARNING! "+self.GetPath()+": Basic datatypes should be completely uppercase.\n"
            if (len(self.Modifiers) > 3): ## a datatype can at most have 3 modifiers
                ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+": A datatype can at most have 3 modifieres (UNIQUE xxx XREF yyy REPEAT).\n"
            else:
                XREF_count=0
                UNIQUE_count=0
                REPEAT_count=0
                for Modif in self.Modifiers:
                    if (isinstance(Modif,UNIQUE_Modifier)):
                        UNIQUE_count=UNIQUE_count+1
                    else:
                        if (isinstance(Modif,XREF_Modifier)):
                            ERR=Modif.Validate(BasicDataTypes)
                            if (ERR != ""):
                                ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+" "+Modif.__str__()+": "+ERR
                            
                            if (self.GetName()[0] != "?"): ## only classes can have references
                                ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+": Only classes can have XREF modifier.\n"
                            if (ClassDict.has_key(self.GetName())):
                                CLASS=ClassDict[self.GetName()]
                                if (CLASS.HasLabel(Modif.GetRef()) == ""): ## checks wether the reference exists
                                    ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+": The label '"+Modif.GetRef()+"' is not defined in class " +self.GetName()+".\n"
                                else:
                                    LABELNAME=string.split(CLASS.HasLabel(Modif.GetRef()))[-1]
                                    if (LABELNAME != Modif.GetRef()):
                                        ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": WARNING! Possible spelling problem of reference '"+Modif.GetRef()+"'"+", probably should be '"+string.split(CLASS.HasLabel(Modif.GetRef()))[-1]+"'.\n"
                                    
                            XREF_count=XREF_count+1
                        else:
                            if (isinstance(Modif,REPEAT_Modifier)):
                                REPEAT_count=REPEAT_count+1
                                if (self.IsLast() !=1): ## only the last datatype can have modifier REPEAT
                                    ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+": Only the last datatype in a line can have a REPEAT modifier.\n"
                if ((REPEAT_count > 1) or (XREF_count > 1) or (UNIQUE_count > 1)): ## each modifier can only appear once
                    ERROR=ERROR+"line "+str(self.GetLineNumber(ClassDict))+": ERROR in "+self.GetPath()+": Each modifier can only be used once for a datatype.\n"
        return ERROR
                

    def AddModifier(self,Modif):      ## adds a modifier
        self.Modifiers.append(Modif)

    def RemoveModifier(self,ModName):    ## removes a modifier
        for Modif in self.Modifiers:
            if (Modif.GetName() == ModName):
                self.Modifiers.remove(Modif)
                

    def SetName(self,Name):       ## sets the name
        self.Name=Name

    def GetName(self):            ## returns the name
        return self.Name

    def GetPath(self):             ## returns the path to the object
        PATH=""
        if (hasattr(self,"Parent")):
            PATH=self.Parent.GetPath()
        if (PATH==""):
            PATH=self.Name
        else:
            PATH=PATH+" -> "+self.Name
        return PATH


    def __str__(self):                 ## returns a string representation of the object
        STRING=self.GetName()
        for Modif in self.Modifiers:
            if (isinstance(Modif,UNIQUE_Modifier)):
                STRING=Modif.__str__()+" "+STRING
            else:
                STRING=STRING+" "+Modif.__str__()
        if (self.IsLast()==1):
            STRING=STRING+"\n"
        return STRING

        

    def __repr__(self):  ## similar to __str__
        return self

    def IsLast(self):     ## checks wether the datatype is the last one in a line
        RESULT=0
        POSITION=self.Parent.FindElement(self)
        if ((POSITION[1]+1) == self.Parent.GetLenPart(POSITION[0])):
            RESULT=1
        return RESULT
    

    def HasLabel(self,LabelName):  ## this is needed to ensure that the recursive method for class and labels works
        return ""

    def GetLineNumber(self,ClassDict):
        MODEL=ClassDict[string.split(self.GetPath())[0]].GetParent()    
        if (isinstance(MODEL,Model)):
            LINENUMBER=MODEL.FindLine(self.GetPath())
        else:
            LINENUMBER=-1
        return LINENUMBER

        
class Modifier:  ## abstract modifier class
    def __init__(self,Parent,Name):
        self.SetParent(Parent)
        self.SetName(Name)

    def __str__(self):
        return self.GetName()

    def __repr__(self):
        return self

    def SetName(self,Name):
        self.Name=Name

    def SetParent(self,Parent):
        self.Parent=Parent

    def GetName(self):
        return self.Name

    def GetParent(self):
        return self.Parent

    def Validate(self,BasicDataTypes):
        return ""
    


class XREF_Modifier(Modifier):    ## inherits from Modifier and adds some extra functinality
    def __init__(self,Parent,Name,Reference):
        Modifier.__init__(self,Parent,Name)
        self.SetRef(Reference)

    def __str__(self):
        return Modifier.__str__(self)+" "+self.GetRef()

    def SetRef(self,Ref):
        self.Reference=Ref

    def GetRef(self):
        return self.Reference

    def Validate(self,BasicDataTypes):
        Error = ""
        if ((self.GetRef()[0] in ["#","?"]) or (string.upper(self.GetRef()) in map(string.upper,BasicDataTypes))):
            Error="ERROR. A reference can not be a class or datatype, it must be a label.\n"
        return Error


class REPEAT_Modifier(Modifier):  ## same as Modifier just needed for convienence
    pass
        

class UNIQUE_Modifier(Modifier):  ## same as Modifier
    pass






if __name__=="__main__":  ## this method is the main method
                          ## it is executed when the file is called from the command line
                          ## it is not executed when the file is imported as a module
    if (len(sys.argv) > 2):
        print "usage ParseModel.py 'modelfile'"

    else:
        File=sys.argv[1]            ## first command line argument is the filename
        Text=open(File,'r').readlines()   ## read all the text from the file
        M=Model(Text)                   ## creates an instance M of the Model
        print M.Validate()              ## prints the result of the validation process for the whole model




