#! /usr/bin/env python

import clang.cindex
import sys
import clang
from clang.cindex import CursorKind, TokenKind, StorageClass

clang.cindex.Config.set_library_file("/usr/lib/llvm-5.0/lib/libclang.so.1")
human_readable=False

def fully_qualified(c):
    if c is None:
        return ''
    elif c.kind == CursorKind.TRANSLATION_UNIT:
        return ''
    else:
        res = fully_qualified(c.semantic_parent)
        if res != '':
            return res + '::' + c.spelling
    return c.spelling

def output(name, start, end):
    if human_readable:
        print "{} {} {}".format(name, start, end)
    else:
        print "{},{}p".format(start,end)


def recurse_into_static_functions(cursor, static_functions_dict):
    for d in cursor.walk_preorder():
        if d.kind == CursorKind.CALL_EXPR and d.referenced.storage_class == StorageClass.STATIC:
                                #print "\t->{}".format(fully_qualified(d.referenced))
            if d.referenced.extent.start.line not in static_functions_dict:
                static_functions_dict[d.referenced.extent.start.line] = d.referenced
                recurse_into_static_functions(d.referenced, static_functions_dict)


if __name__ == "__main__":
    idx = clang.cindex.Index.create()
    tu = idx.parse(sys.argv[1])

    function_dict = {}
    kindtypes = dict()
    static_functions_dict = dict()
    for idx, c in enumerate(tu.cursor.walk_preorder()):

        kindtypes[c.kind] = kindtypes.get(c.kind, 0) + 1

        if c.kind == CursorKind.FUNCTION_DECL and c.get_definition() is not None:
            #print("\n".join(dir(c)))
            #sys.exit(0)
            if str(c.location.file).endswith("c"):

                if c.storage_class == StorageClass.STATIC:
                    #print "{} : {}--{}".format(fully_qualified(c.referenced), c.extent.start.line,
                    #
                    #c.extent.end.line)
                    pass

                with open(sys.argv[2], "r") as headerfile:
                    if "{}(".format(c.referenced.spelling) in headerfile.read():
                        #print "{} : {}--{}".format(fully_qualified(c.referenced), c.extent.start.line, c.extent.end.line)

                        function_dict[c.extent.start.line] = c

                        recurse_into_static_functions(c, static_functions_dict)



    #print "\n".join(map(str, kindtypes.items()))
    #
    # loop over documentation that and keep merge tokens that are doxygen comments directly before
    # the start of a new function declaration
    #
    for c in tu.cursor.get_tokens():

        # stop at every comment token
        if c.kind == TokenKind.COMMENT:
            nextlinenumber = c.extent.end.line + 1

            #
            # is this function part of the current header module?
            #
            if  nextlinenumber in function_dict or nextlinenumber in static_functions_dict:

                # uncomment to print the documentation
                #print c.spelling
                thedict = function_dict if nextlinenumber in function_dict else static_functions_dict
                function = thedict[nextlinenumber]
                functionname = fully_qualified(function.referenced)

                #
                # print the whole extent from the beginning to the end.
                #
                output(functionname, c.extent.start.line, function.extent.end.line)



