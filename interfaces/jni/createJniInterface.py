#!/usr/bin/env python

import sys
import os.path
from  os import system

import xml.dom.minidom
from xml.dom.minidom import Node

import string
import pprint
import textwrap

def createText(npara, paratext, kind):

    if npara == 0:
        if kind == "":
            text = textwrap.fill(paratext, initial_indent='   /** ', subsequent_indent='    *  ', width=120) + "\n"
        else:
            initial = '   /** @'  + kind + ' '
            text = textwrap.fill(paratext, initial_indent=initial, subsequent_indent='    *  ', width=120) + "\n"
    else:
        text = "    *\n"
        if kind == "":
            text += textwrap.fill(paratext, initial_indent='    *  ', subsequent_indent='    *  ', width=120) + "\n"
        else:
            width = len(kind)
            initial = '    *  @'  + kind + ' '
            subsequent = '    *'.ljust(width + 9, ' ')
            text += textwrap.fill(paratext, initial_indent=initial, subsequent_indent=subsequent, width=120) + "\n"

    return text


def getText(nodelist):
    rc = []
    for node in nodelist:
        if node.nodeType == node.TEXT_NODE:
            rc.append(node.data)
        else:
            kind = node.getAttribute('kind')
            if kind == "":
                rc.append(getText(node.childNodes))

    text = ''.join(rc)

    #remove leading and tailing white space
    text = text.strip();

    return text

def getNodeText(node):
    rc = []

    if node.nodeType == node.TEXT_NODE:
        rc.append(node.data)
    else:
        rc.append(getText(node.childNodes))

    text = ''.join(rc)

    #remove leading and tailing white space
    text = text.strip();

    return text

# parameter which currently cannot be handled; each method which contains such an parameter gets ignored
blacklistparams = ['SCIP_QUADELEM', 'SCIP_NLPSTATISTICS', 'SCIP_BIVAR_CONVEXITY', 'SCIP_INTERVAL', 'SCIP_EXPR', 'va_list']

# list of methods where the one reference gets not moved to the return value since it used for freeing an object
freemthodlist = []

def getTaglist(node, tag):
    taglist = node.getElementsByTagName(tag)
    return taglist

def getTextFromTag(node, tag):
    taglist = node.getElementsByTagName(tag)

    if taglist.length > 0:
        textnode = taglist[0];
        texttag = getText(textnode.childNodes)
    else:
        texttag = ""

    return texttag

# parse parameters with name description
def parseParameters(param):

    params = []

    for param in param.getElementsByTagName('param'):

        # parameter name
        paramlist = param.getElementsByTagName('declname')
        if paramlist.length > 0:
            textnode = paramlist[0]
            paramname = "j" + getText(textnode.childNodes)
        else:
            paramname = ""

        # parameter type
        paramlist = param.getElementsByTagName('type')
        if paramlist.length > 0:
            textnode = paramlist[0]
            paramtype = getText(textnode.childNodes)
        else:
            paramtype = ""

        # parameter description
        paramdesc = getTextFromTag(param, 'briefdescription')

        params.append((paramname, paramtype, paramdesc))

    return params

# get file class name
def getComponent(doc):
    compoundlist = doc.getElementsByTagName('compoundname')

    if compoundlist.length > 0:
        filename = getText(compoundlist[0].childNodes)

        filename = filename.replace('.h','')

        if filename.startswith("pub_"):
            filename = filename.replace('pub_','').capitalize()
            return filename

        if filename.count('_') == 1:
            namelist = filename.split('_')
            namelist = [name.capitalize() for name in namelist]
            filename = ''.join(namelist)
            return filename

        if filename.startswith("scip"):
            return ""

    return "-1"

# basis data types
def convertBasicType(basictype):
    basictype = basictype.replace('SCIP_Bool', 'boolean')
    basictype = basictype.replace('SCIP_Real', 'double')
    basictype = basictype.replace('SCIP_Longint', 'long')
    basictype = basictype.replace('const', '')
    basictype = basictype.replace('char *', 'String')
    basictype = basictype.replace('unsigned int', 'int')
    return basictype

# converts enums pointer to long
def convertEnumsType(enumtype):
    # grep "typedef enum" src/*/*h | sed 's/[ ]\+/ /g' | cut -d " " -f 4 | sed "s/;/',/" | sed '/^$/d' | sed "s/^/'/"
    # test = subprocess.check_output(["grep \"typedef enum\" src/*/*h | sed 's/[ ]\+/ /g' | cut -d \" \" -f 4 | sed \"s/;/',/\" | sed '/^$/d' | sed \"s/^/'/\""])

    enumtypes = [
        'SCIP_EXPROP',
        'SCIP_EXPRCURV',
        'SCIP_NLPPARAM',
        'SCIP_NLPSOLSTAT',
        'SCIP_NLPTERMSTAT',
        'SCIP_LINEARCONSTYPE',
        'SCIP_SETPPCTYPE',
        'SCIP_CLOCKTYPE',
        'SCIP_DISPSTATUS',
        'SCIP_BRANCHDIR',
        'SCIP_LPSOLSTAT',
        'SCIP_BOUNDTYPE',
        'SCIP_SIDETYPE',
        'SCIP_ROWORIGINTYPE',
        'SCIP_LPALGO',
        'SCIP_OBJSEN',
        'SCIP_LPPARAM',
        'SCIP_PRICING',
        'SCIP_BASESTAT',
        'SCIP_LPSOLQUALITY',
        'SCIP_VERBLEVEL',
        'SCIP_PARAMTYPE',
        'SCIP_PARAMSETTING',
        'SCIP_PARAMEMPHASIS',
        'SCIP_OBJSENSE',
        'SCIP_RESULT',
        'SCIP_RETCODE',
        'SCIP_EFFICIACYCHOICE',
        'SCIP_STAGE',
        'SCIP_SETTING',
        'SCIP_SOLORIGIN',
        'SCIP_STATUS',
        'SCIP_NODETYPE',
        'SCIP_VARSTATUS',
        'SCIP_VARTYPE',
        'SCIP_DOMCHGTYPE',
        'SCIP_BOUNDCHGTYPE',
        'SCIP_VBCCOLOR',
        'TCLIQUE_STATUS',
        'SCIP_LINEARCONSTYPE'
        ]

    for enum in reversed(enumtypes):
        enumtype = enumtype.replace(enum, 'int')
    return enumtype

# converts bit arrays to int
def convertUnsignedInt(unsignedinttype):
    # grep "typedef unsigned int" src/*/*h | sed 's/[ ]\+/ /g' | cut -d " " -f 4 | sed "s/;/',/" | sed '/^$/d' | sed "s/^/'/"
    unsignedinttypes = [
        'SCIP_EXPRINTCAPABILITY',
        'SCIP_SINGLEPACKET',
        'SCIP_DUALPACKET',
        'SCIP_EVENTTYPE',
        'SCIP_PROPTIMING',
        'SCIP_PRESOLTIMING',
        'SCIP_HEURTIMING'
        ]
    for unsignedint in unsignedinttypes:
        unsignedinttype = unsignedinttype.replace(unsignedint, 'int')

    return unsignedinttype

# converts struct pointer to long
def convertStructPointer(pointertype):
    # grep "typedef struct" src/*/*h | sed 's/[ ]\+/ /g' | cut -d " " -f 4 | sed "s/;/\ \*',/" | sed '/^$/d' | sed "s/^/'/"
    # grep "typedef union" src/*/*h | sed 's/[ ]\+/ /g' | cut -d " " -f 4 | sed "s/;/\ \*',/" | sed '/^$/d' | sed "s/^/'/"
    pointertypes =[
        'BMS_CHKMEM *',
        'BMS_BLKMEM *',
        'DIJKSTRA_GRAPH *',
        'SCIP_NLPIORACLE *',
        'SCIP_EXPR *',
        'SCIP_EXPRTREE *',
        'SCIP_QUADELEM *',
        'SCIP_EXPRDATA_QUADRATIC *',
        'SCIP_EXPRDATA_MONOMIAL *',
        'SCIP_EXPRDATA_POLYNOMIAL *',
        'SCIP_EXPRGRAPHNODE *',
        'SCIP_EXPRGRAPH *',
        'SCIP_EXPRINT *',
        'SCIP_EXPRINTDATA *',
        'SCIP_NLPI *',
        'SCIP_NLPIDATA *',
        'SCIP_NLPIPROBLEM *',
        'SCIP_NLPSTATISTICS *',
        'SCIP_LINCONSUPGRADE *',
        'SCIP_QUADVAREVENTDATA *',
        'SCIP_QUADVARTERM *',
        'SCIP_BILINTERM *',
        'SCIP_INTERVAL *',
        'SCIP_FILE *',
        'SCIP_BOOLPARAM *',
        'SCIP_INTPARAM *',
        'SCIP_LONGINTPARAM *',
        'SCIP_REALPARAM *',
        'SCIP_CHARPARAM *',
        'SCIP_STRINGPARAM *',
        'SCIP_BRANCHCAND *',
        'SCIP_BRANCHRULE *',
        'SCIP_BRANCHRULEDATA *',
        'SCIP_BUFFER *',
        'SCIP_CLOCK *',
        'SCIP_CPUCLOCK *',
        'SCIP_WALLCLOCK *',
        'SCIP_CONFLICTHDLR *',
        'SCIP_CONFLICTHDLRDATA *',
        'SCIP_CONFLICTSET *',
        'SCIP_LPBDCHGS *',
        'SCIP_CONFLICT *',
        'SCIP_CONSHDLR *',
        'SCIP_CONS *',
        'SCIP_CONSHDLRDATA *',
        'SCIP_CONSDATA *',
        'SCIP_CONSSETCHG *',
        'SCIP_CUTPOOL *',
        'SCIP_CUT *',
        'SCIP_DIALOG *',
        'SCIP_DIALOGDATA *',
        'SCIP_DIALOGHDLR *',
        'SCIP_LINELIST *',
        'SCIP_DISP *',
        'SCIP_DISPDATA *',
        'SCIP_EVENTHDLR *',
        'SCIP_EVENTHDLRDATA *',
        'SCIP_EVENT *',
        'SCIP_EVENTVARADDED *',
        'SCIP_EVENTVARDELETED *',
        'SCIP_EVENTVARFIXED *',
        'SCIP_EVENTVARUNLOCKED *',
        'SCIP_EVENTOBJCHG *',
        'SCIP_EVENTBDCHG *',
        'SCIP_EVENTHOLE *',
        'SCIP_EVENTIMPLADD *',
        'SCIP_EVENTROWADDEDSEPA *',
        'SCIP_EVENTROWDELETEDSEPA *',
        'SCIP_EVENTROWADDEDLP *',
        'SCIP_EVENTROWDELETEDLP *',
        'SCIP_EVENTROWCOEFCHANGED *',
        'SCIP_EVENTROWCONSTCHANGED *',
        'SCIP_EVENTROWSIDECHANGED *',
        'SCIP_EVENTDATA *',
        'SCIP_EVENTFILTER *',
        'SCIP_EVENTQUEUE *',
        'SCIP_HEUR *',
        'SCIP_HEURDATA *',
        'SCIP_HISTORY *',
        'SCIP_VALUEHISTORY *',
        'SCIP_VBOUNDS *',
        'SCIP_IMPLICS *',
        'SCIP_CLIQUE *',
        'SCIP_CLIQUETABLE *',
        'SCIP_CLIQUELIST *',
        'SCIP_INTERRUPT *',
        'SCIP_COLSOLVALS *',
        'SCIP_ROWSOLVALS *',
        'SCIP_LPSOLVALS *',
        'SCIP_COL *',
        'SCIP_ROW *',
        'SCIP_LP *',
        'SCIP_LPI *',
        'SCIP_LPISTATE *',
        'SCIP_MEM *',
        'SCIP_MESSAGEHDLR *',
        'SCIP_MESSAGEHDLRDATA *',
        'SCIP_SPARSESOL *',
        'SCIP_QUEUE *',
        'SCIP_PQUEUE *',
        'SCIP_QUEUE *',
        'SCIP_HASHTABLE *',
        'SCIP_HASHTABLELIST *',
        'SCIP_HASHMAP *',
        'SCIP_HASHMAPLIST *',
        'SCIP_REALARRAY *',
        'SCIP_INTARRAY *',
        'SCIP_BOOLARRAY *',
        'SCIP_PTRARRAY *',
        'SCIP_RESOURCEACTIVITY *',
        'SCIP_PROFILE *',
        'SCIP_RESOURCEACTIVITY *',
        'SCIP_DIGRAPH *',
        'SCIP_BT *',
        'SCIP_BTNODE *',
        'SCIP_NLROW *',
        'SCIP_NLP *',
        'SCIP_NODEPQ *',
        'SCIP_NODESEL *',
        'SCIP_NODESELDATA *',
        'SCIP_PARAM *',
        'SCIP_PARAMDATA *',
        'SCIP_PARAMSET *',
        'SCIP_PRESOL *',
        'SCIP_PRESOLDATA *',
        'SCIP_PRICER *',
        'SCIP_PRICERDATA *',
        'SCIP_PRICESTORE *',
        'SCIP_PRIMAL *',
        'SCIP_PROB *',
        'SCIP_PROBDATA *',
        'SCIP_PROP *',
        'SCIP_PROPDATA *',
        'SCIP_READER *',
        'SCIP_READERDATA *',
        'SCIP_RELAX *',
        'SCIP_RELAXATION *',
        'SCIP_RELAXDATA *',
        'SCIP *',
        'SCIP_SEPA *',
        'SCIP_SEPADATA *',
        'SCIP_SEPASTORE *',
        'SCIP_SET *',
        'SCIP_SOL *',
        'SCIP_STAT *',
        'SCIP_PROBINGNODE *',
        'SCIP_SIBLING *',
        'SCIP_CHILD *',
        'SCIP_LEAF *',
        'SCIP_JUNCTION *',
        'SCIP_PSEUDOFORK *',
        'SCIP_FORK *',
        'SCIP_SUBROOT *',
        'SCIP_NODE *',
        'SCIP_PENDINGBDCHG *',
        'SCIP_TREE *',
        'SCIP_DOMCHGBOUND *',
        'SCIP_DOMCHGBOTH *',
        'SCIP_DOMCHGDYN *',
        'SCIP_BOUNDCHG *',
        'SCIP_BDCHGIDX *',
        'SCIP_BDCHGINFO *',
        'SCIP_BRANCHINGDATA *',
        'SCIP_INFERENCEDATA *',
        'SCIP_HOLECHG *',
        'SCIP_HOLE *',
        'SCIP_HOLELIST *',
        'SCIP_DOM *',
        'SCIP_ORIGINAL *',
        'SCIP_AGGREGATE *',
        'SCIP_MULTAGGR *',
        'SCIP_NEGATE *',
        'SCIP_VAR *',
        'SCIP_VARDATA *',
        'SCIP_VBC *',
        'TCLIQUE_GRAPH *',
        'TCLIQUE_DATA *',
        'XML_ATTR *',
        'XML_NODE *',
        'SCIP_EXPROPDATA *',
        'SCIP_DOMCHG *',
        'void *',
        'FILE *'
        ]
    for pointer in pointertypes:
        pointertype = pointertype.replace(pointer, 'long', 1)
        
    return pointertype

def convertArray(arraytype):
    arraytype = arraytype.replace('*', '[]')
    return arraytype

def convertFunction(functiontype):
    return functiontype

def convertReturnType(returntype):
    returntype = convertBasicType(returntype)
    returntype = convertEnumsType(returntype)
    returntype = convertStructPointer(returntype)
    returntype = convertUnsignedInt(returntype)
    returntype = convertArray(returntype)
    return returntype

def convertParamType(paramtype):
    paramtype = convertBasicType(paramtype)
    paramtype = convertEnumsType(paramtype)
    paramtype = convertStructPointer(paramtype)
    paramtype = convertUnsignedInt(paramtype)
    paramtype = convertArray(paramtype)
    paramtype = convertFunction(paramtype)
    return paramtype

def createEnums(filename):
    print "create Enum" + filename

    if os.path.exists(filename) and os.path.isfile(filename):
        filename = os.path.abspath(filename)

    # parse the xml file
    doc = xml.dom.minidom.parse(filename)

    # loop over all methods and add them to the interface
    for method in doc.getElementsByTagName("memberdef"):

        # docu
        kind = method.getAttribute("kind")
        if kind != 'enum':
            continue

        enumname = getTextFromTag(method, 'name')

        if enumname.startswith("SCIP_"):
            enumname = enumname.replace('SCIP_','').capitalize()

        if filename.count('_') == 1:
            namelist = enumname.split('_')
            namelist = [name.capitalize() for name in namelist]
            enumname = ''.join(namelist)

        print "enum " + enumname

        jnifile = open("java/de/zib/jscip/nativ/jni/JniScip" + enumname + ".java", 'w')

        jnifile.write("/* DO NOT EDIT THIS FILE - it is machine generated */\n")
        jnifile.write("\n")
        jnifile.write("package de.zib.jscip.nativ.jni;\n")
        jnifile.write("\n")
        jnifile.write("public class JniScip" + enumname  + " {\n")
        jnifile.write("\n")

        # loop over all enum values
        for enum in method.getElementsByTagName("enumvalue"):

            enumname = getTextFromTag(enum, 'name')
            initializer = getTextFromTag(enum, 'initializer')
            enumdetaileddescriptions = []

            taglist = getTaglist(enum, 'detaileddescription')

            if taglist.length > 1:
                print "Error"
                exit

            detaileddescriptions = []

            # parse detailed descriptions
            if taglist.length > 0:

                npara = 0
                paralist = getTaglist(taglist[0], 'para')

                # loop over all paragraphs
                for para in taglist[0].childNodes:
                    if para.nodeType == para.TEXT_NODE:
                        continue

                    simplesectlsit = getTaglist(para, 'simplesect')

                    for simplesect in simplesectlsit:
                        kind = simplesect.getAttribute('kind')
                        paratext = getNodeText(simplesect);
                        detaileddescriptions.append(createText(npara, paratext, kind))
                        npara += 1

                    paratext = getNodeText(para);
                    if paratext != "":
                        detaileddescriptions.append(createText(npara, paratext, ""))
                        npara += 1

            # print description
            for detaileddescription in detaileddescriptions:
                jnifile.write(detaileddescription)

            jnifile.write("    */\n")
            jnifile.write("   final public static int " + enumname + " = " + initializer.rsplit(' ')[-1] + ";\n");
            jnifile.write("\n")

        jnifile.write("}\n")



#################################################
#
# Main method
#
#################################################

# start if the main loop
if len(sys.argv) < 2:
    print "Usage: %s [<xml-doxygen-filename>, ...]" % sys.argv[0]
    sys.exit(1)

# remove NativeScip.java and JniScip.java
if os.path.isfile("java/de/zib/jscip/nativ/NativeScip.java") :
    os.remove("java/de/zib/jscip/nativ/NativeScip.java");

if os.path.isfile("java/de/zib/jscip/nativ/jni/JniScip.java") :
    os.remove("java/de/zib/jscip/nativ/jni/JniScip.java");

# create for each imput file one output file
for arg in sys.argv[1:]:

    # check if th file name contains "8h" which indicates an header file
    if not "8h" in arg:
        continue

    # ignore all type headers for now
    if "type" in arg:
        createEnums(arg)
        continue

    # ignore all C++ wrapper classes
    if "/obj" in arg:
        continue

    # ignore all file io methods
    if "fileio" in arg:
        continue

    print "*****************************************"
    print "* " + arg
    print "*****************************************"

    if os.path.exists(arg) and os.path.isfile(arg):
        arg = os.path.abspath(arg)

    # parse the xml file
    doc = xml.dom.minidom.parse(arg)

    filename = getTextFromTag(doc, 'compoundname')
    component = getComponent(doc)

    # unknown file skip that file completely
    if component == "-1":
        print "WARNING: unknown file"
        continue

    # if there is no .c implementation, we should skip the complete file
    if not os.path.isfile("src/JniScip" + component + ".c"):
        print "File src/JniScip" + component + ".c, does not exist. Skip it"
        continue

    # check if the java files already exist
    if component == "" and os.path.isfile("java/de/zib/jscip/nativ/NativeScip.java") :
        nativejnifile = open("java/de/zib/jscip/nativ/NativeScip.java")
        nativelines = nativejnifile.readlines()
        nativejnifile.close()

        jnifile = open("java/de/zib/jscip/nativ/jni/JniScip.java")
        jnilines = jnifile.readlines()
        jnifile.close()
    else:
        nativelines = []
        jnilines = []

    if not os.path.exists("java/de/zib/jscip/nativ/jni"):
        os.makedirs("java/de/zib/jscip/nativ/jni")

    nativejnifile = open("java/de/zib/jscip/nativ/NativeScip" + component + ".java", 'w')
    jnifile = open("java/de/zib/jscip/nativ/jni/JniScip" + component + ".java", 'w')
    cfile = open("src/JniScip" + component + ".c", 'r')

    # create preamble for Native interface
    if nativelines == []:
        nativejnifile.write("/* DO NOT EDIT THIS FILE - it is machine generated */\n")
        nativejnifile.write("\n")
        nativejnifile.write("package de.zib.jscip.nativ;\n")
        nativejnifile.write("\n")
        nativejnifile.write("/** JNI interface for " + filename + " */\n")
        nativejnifile.write("public interface NativeScip" + component + " {\n")
        nativejnifile.write("\n")
    else:
        nativejnifile.writelines([item for item in nativelines[:-1]])


    # create preamble for JNI (native) implementation
    if jnilines == []:
        jnifile.write("/* DO NOT EDIT THIS FILE - it is machine generated */\n")
        jnifile.write("\n")
        jnifile.write("package de.zib.jscip.nativ.jni;\n")
        jnifile.write("\n")
        jnifile.write("import de.zib.jscip.nativ.NativeScip" + component +";\n")
        jnifile.write("import de.zib.jscip.nativ.NativeScipException;\n")
        jnifile.write("\n")
        jnifile.write("/** JNI implementation (native) for " + filename + " */\n")
        jnifile.write("public class JniScip" + component + " implements NativeScip" + component + "{\n")
        jnifile.write("\n")
    else:
        jnifile.writelines([item for item in jnilines[:-1]])

    # loop over all methods and add them to the interface
    for method in doc.getElementsByTagName("memberdef"):

        # docu
        kind = method.getAttribute("kind")
        if kind != 'function':
            continue

#        if method.hasAttributes():
#            attributes = method.attributes()
#            print attributes

        methodname = getTextFromTag(method, 'name')

        # if a function is called  fun_name, javah changes the name to fun_1name
        javamethodname = "(" + methodname.replace('SCIP', '').replace('_','_1') + ")"
        cfile.seek(0)
        skipmethod = 1
        for line in cfile:
            if javamethodname in line:
                skipmethod = 0
                break
        if skipmethod:
            print "Skip method " + methodname + " not implemented on " + cfile.name + " " + javamethodname + " not found"
            continue;

        # check if the method which potentially frees an object
        if "free" in methodname or "release" in methodname:
            freemethod = 1;
        else:
            freemethod = 0;

        taglist = getTaglist(method, 'detaileddescription')

        if taglist.length > 1:
            print "Error"
            exit

        detaileddescriptions = []

        # parse detailed descriptions
        if taglist.length > 0:

            npara = 0
            paralist = getTaglist(taglist[0], 'para')

            # loop over all paragraphs
            for para in taglist[0].childNodes:
                if para.nodeType == para.TEXT_NODE:
                    continue

                simplesectlsit = getTaglist(para, 'simplesect')

                for simplesect in simplesectlsit:
                    kind = simplesect.getAttribute('kind')
                    paratext = getNodeText(simplesect);
                    detaileddescriptions.append(createText(npara, paratext, kind))
                    npara += 1

                paratext = getNodeText(para);
                if paratext != "":
                    detaileddescriptions.append(createText(npara, paratext, ""))
                    npara += 1

        # parse return type
        returntype = getTextFromTag(method, 'type')

        # method
        nreturnparams = 0
        nfunctionparam = 0
        skipmethod = 0

        # check if the return type is black listed; if this is the case ignore the method completely
        for  blacklistparam in blacklistparams:
            if returntype.find(blacklistparam) != -1:
                skipmethod = 1
                break;

        # parse parameters
        params = parseParameters(method)

        # loop of the parameters and check them for references, create the parameter list as string, ...
        for param in params:
            if skipmethod == 1:
                break;

            # currently we cannot handle the "..." argument
            if param[1] == "...":
                skipmethod = 1
                break

            # check if the parameter is black listed; if this is the case ignore the method completely
            for  blacklistparam in blacklistparams:
                if param[1].find(blacklistparam) != -1:
                    skipmethod = 1
                    break;

            # check if the parameters is return parameter (reference)
            if param[2].startswith('pointer to') or param[2].startswith('stores') or param[2].startswith('array to store'):
                nreturnparams += 1
                continue

            # check if the parameters is a function pointer
            if param[1].startswith('SCIP_DECL_'):
                nfunctionparam += 1
                continue

            if param[1] == "void":
                continue

        # in case we have only one reference within the parameter list, the
        # return type is an SCIP_RETCODE, and the method is not a "free"
        # method, we move that reference to the return type;
        if nreturnparams == 1 and returntype == "SCIP_RETCODE" and freemethod == 0:
            # remove that one return parameter from the parameter list and use it as return value
            for param in params:
                if param[2].startswith('pointer to') or param[2].startswith('stores'):

                    # remove the one of the stars
                    returntype = param[1].replace('*', '', 1)
                    returndesc = param[2].replace('pointer to', '')
                    returndesc = param[2].replace('stores', '')
                    nreturnparams = 0
                    params.remove(param)
                    break

                if param[2].startswith('array to store'):

                    # remove the one of the stars
                    returntype = param[1].replace('*', '', 1) + '[]'
                    print returntype
                    returndesc = param[2].replace('array to store', '')
                    nreturnparams = 0
                    params.remove(param)
                    break


        if nreturnparams == 1 and freemethod == 1:
            for param in params:
                if param[2].startswith('pointer to') or param[2].startswith('stores'):
                    # remove the one of the stars
                    paramname = param[0]
                    paramtype = param[1].replace('*', '', 1)
                    paramdesc = param[2].replace('pointer to', '')
                    params.remove(param)
                    params.append((paramname, paramtype, paramdesc))
                    nreturnparams = 0
                    break

        if skipmethod == 1:
            continue

        if nfunctionparam > 0:
            print "ignoring " + methodname + " since it contains a function pointer"
            continue

        if nreturnparams > 0:
            print "ignoring " + methodname + " since it contains more than one return value"
            continue

        if returntype == "":
            print "ignoring " + methodname + " due to unknown return type"
            continue

        javaparamstr = []
        njavaparams = 0

        # create parameter string
        for param in params:
            if param[1] == "void":
                continue

            if njavaparams > 0:
                javaparamstr.append(', ')

            paramtype = convertParamType(param[1])
            javaparamstr.append(paramtype + " " + param[0])
            njavaparams += 1

        # print description
        for detaileddescription in detaileddescriptions:
            nativejnifile.write(detaileddescription)
            jnifile.write(detaileddescription)

        # print parameter list
        paramformat="    *  @param %-20s   %s\n"

        nparams = 0

        for param in params:
            nparams += 1

            if param[1] != "void":
                if nparams == 1:
                    nativejnifile.write("    *\n")
                    jnifile.write("    *\n")

                nativejnifile.write(paramformat % (param[0], param[2]))
                jnifile.write(paramformat % (param[0], param[2]))

        methodname = methodname.replace('SCIP', '')

        if returntype == "SCIP_RETCODE":
            returntype = "void"

        returntype = convertReturnType(returntype)

        nativejnifile.write("    */\n")
        nativejnifile.write("   public " + returntype + " " + methodname + "(" + ''.join(javaparamstr) + ")\n")
        nativejnifile.write("      throws NativeScipException;\n")

        jnifile.write("    */\n")
        jnifile.write("   public native " + returntype + " " + methodname + "(" + ''.join(javaparamstr) + ")\n")
        jnifile.write("      throws NativeScipException;\n")

        nativejnifile.write("\n")
        jnifile.write("\n")

    nativejnifile.write("}\n")
    jnifile.write("}\n")

    nativejnifile.close()
    jnifile.close()
