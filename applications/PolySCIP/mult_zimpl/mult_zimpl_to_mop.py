#!/usr/bin/python3

import argparse
from collections import defaultdict
import os
import re
import subprocess
import sys

def obj_name_is_feasible(obj, names):
    if len(obj) > 6:
        print("Name of objective is too long. Please make sure " \
              "that all names of the objectives are less than 7 characters.")
        return False
    elif obj in names:
        print("Duplicate names of objectives. Please make sure that " \
              "all names of the objectives are unique.")
        return False

    return True


def zimpl_to_mop(problem_file, path_for_mop_file, output_basename, path_to_zimpl):
    """Takes a multicriteria problem file in extended zimpl format and
    adjusts it such that it conforms entirely to the single-criteria zimpl format. 
    Zimpl is then called on the file and its output adjusted to the mop file format.
    """
    min_max_pattern = re.compile(r'[ \t]*m(ini|axi)mize[ \t]+[a-z]\w*[ \t]*:', re.IGNORECASE) 
    obj_pattern = re.compile(r'[\t]*[a-z]\w*[ \t]*:', re.IGNORECASE)
    subto_pattern = re.compile(r'[ \t]*subto[ \t]*:')

    if path_for_mop_file and not path_for_mop_file.endswith("/"):
        path_for_mop_file += "/"
    if path_to_zimpl and not path_to_zimpl.endswith("/"):
        path_to_zimpl += "/"

    problem_filename = os.path.basename(problem_file)
    splitted = problem_filename.split('.') 
    basename, suffix = splitted[0], splitted[-1]
    file_suffix = ".mzpl"
    if file_suffix == "."+suffix:
        file_suffix = ".multzpl"

    mzpl = path_for_mop_file
    if output_basename is None:
        mzpl += basename
    else:
        mzpl += output_basename
    objectives = []
    min_max_pattern_found = False

    with open(problem_file, mode='r', encoding='utf-8') as src_file:
        with open(mzpl+file_suffix, mode='w', encoding='utf-8') as mzpl_file:

            for line in src_file:
                if not min_max_pattern_found:
                    min_max_match = min_max_pattern.match(line)
                    if min_max_match:
                        min_max_pattern_found = True
                        substr = min_max_match.group().lstrip()
                        obj_sense = substr[:8].lower()
                        obj_name = substr[8:len(substr)-1].rstrip()
                        if not obj_name_is_feasible(obj_name, objectives):
                            sys.exit()
                        else:
                            objectives.append(obj_name)
                            line = line.replace(obj_sense, "subto", 1)
                            line = line.replace(";", " == 0;", 1)

                else:
                    subto_pattern_match = subto_pattern.match(line)
                    if not subto_pattern_match:
                        obj_pattern_match = obj_pattern.match(line)
                        if obj_pattern_match:
                            substr = obj_pattern_match.group().lstrip()
                            obj_name = substr[:len(substr)-1].rstrip()
                            if not obj_name_is_feasible(obj_name, objectives):
                                sys.exit()
                            else:
                                objectives.append(obj_name)
                                line = "subto " + line
                                line = line.replace(";", " == 0;", 1)

                mzpl_file.write(line)

            proc = subprocess.Popen([path_to_zimpl+"zimpl", '-tmps',
                                     '-o'+mzpl, '-v0', mzpl_file.name])
    os.wait()

    zimpl_obj = "[ \t]*N[ \t]+OBJECTIV"
    zimpl_obj_pattern = re.compile(zimpl_obj)
    comment = "[ \t]*\*"
    comment_pattern = re.compile(comment)
    rows = "[ \t]*ROWS"
    rows_pattern = re.compile(rows)
    if objectives:
        objs = "[ \t]*E[ \t]+(" + "|".join([x for x in objectives]) + ")"
    else:
        objs = r'(?!)'
    objs_pattern = re.compile(objs)


    with open(mzpl+".mps", mode='r', encoding='utf-8') as mps_file:
        with open(mzpl+".mop", mode='w', encoding='utf-8') as mop_file:
            for line in mps_file:
                zimpl_obj_match = zimpl_obj_pattern.match(line)
                comment_match = comment_pattern.match(line)
                row_match = rows_pattern.match(line)
                objs_match = objs_pattern.match(line)
                if zimpl_obj_match or comment_match:
                    continue
                elif row_match:
                    line = "OBJSENSE \n " + obj_sense[0:3].upper() + "\n" + line
                elif objs_match:
                    line = line.replace("E", "N", 1)
                mop_file.write(line)


    os.remove(mzpl+".tbl")
    os.remove(mzpl+file_suffix)
    os.remove(mzpl+".mps")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Process multicriteria zimpl file")
    parser.add_argument("file", help="file containing the multicriteria problem", type=str)
    parser.add_argument("-p", help="path where the output file should be saved; default is ./", metavar="PATH", type=str, default="./")
    parser.add_argument("-o", metavar="BASENAME",
                        help="basename for the output file; default is input file without extension", type=str)                        
    parser.add_argument("--path_to_zimpl", help="path to the zimpl binary", metavar="PATH",
                        type=str, default="")
    args = parser.parse_args()
    zimpl_to_mop(args.file, args.p, args.o, args.path_to_zimpl)
