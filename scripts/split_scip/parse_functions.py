#! /usr/bin/env python

import clang.cindex
import sys
import clang
from clang.cindex import CursorKind, TokenKind, StorageClass

clang.cindex.Config.set_library_file("/usr/lib/llvm-5.0/lib/libclang.so.1")

gap_file_name = "gaps_file.txt"

import argparse

parser = argparse.ArgumentParser(description='Parse input module and header file to create a list of relevant sections')

parser.add_argument('--write_gaps', action='store_true', default=False, dest="write_gaps")
parser.add_argument('files', nargs=2, help='c file and header file')


args = parser.parse_args()

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

def recurse_into_static_functions(cursor, static_functions_dict):
    for d in cursor.walk_preorder():
        if d.kind == CursorKind.CALL_EXPR and d.referenced.storage_class == StorageClass.STATIC:
                                #print "\t->{}".format(fully_qualified(d.referenced))
            if d.referenced.extent.start.line not in static_functions_dict:
                static_functions_dict[d.referenced.extent.start.line] = d.referenced
                recurse_into_static_functions(d.referenced, static_functions_dict)


if __name__ == "__main__":

    header_functions_idx = clang.cindex.Index.create()
    tu  = header_functions_idx.parse(args.files[1], args=['-I../../src'])
    header_functions = set()
    for c in tu.cursor.walk_preorder():
        if c.kind == CursorKind.FUNCTION_DECL and str(c.location.file) == args.files[1]:
            header_functions.add(c.spelling)

    idx = clang.cindex.Index.create()
    tu = idx.parse(args.files[0], args=['-I../../src'])

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


                if c.referenced.spelling in header_functions:
                    #print "{} : {}--{}".format(fully_qualified(c.referenced), c.extent.start.line, c.extent.end.line)

                    function_dict[c.extent.start.line] = c

                    recurse_into_static_functions(c, static_functions_dict)



    #print "\n".join(map(str, kindtypes.items()))
    #
    # loop over documentation that and keep merge tokens that are doxygen comments directly before
    # the start of a new function declaration
    #
    lines = []
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
                #output(functionname, c.extent.start.line, function.extent.end.line)
                lines.append((c.extent.start.line, function.extent.end.line))

    if args.write_gaps:
        with open(gap_file_name, "w") as gap_file:
            start = 1
            i = 0
            gaplist = []
            while i < len(lines) - 1:
                _ , start = lines[i]
                end , _ = lines[i+1]
                if start + 1 <= end - 1:
                    gaplist.append((start + 1, end - 1))
                i += 1

            for (start, end) in gaplist:
                if start <= end:
                    gap_file.write("{} {}\n".format(start, end))
    elif len(lines) > 0:

        #
        # try to open a possibly created gap file
        #
        with open(gap_file_name, "r") as gap_file:
            gaps = {}
            for  line in gap_file:
                gaps_line = map(int, line.split())
                gaps[gaps_line[0]] = gaps_line[1]
            start,end = lines[0]
            end+= 1
            merged_lines = []
            i = 0
            while i < len(lines) - 1:
                nextstart,nextend = lines[i+1]
                nextend += 1
                if end in gaps and gaps[end] == nextstart - 1:
                    end = nextend
                else:
                    merged_lines.append((start,end))
                    start,end = (nextstart,nextend)
                i += 1
            merged_lines.append((start,end - 1))
            lines = merged_lines

        for (start,end) in lines:
            print "{},{}p".format(start,end)
