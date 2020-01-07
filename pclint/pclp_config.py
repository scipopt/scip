#!/usr/bin/python

##############################################################################
# Copyright Gimpel Software LLC 2019. All rights reserved.
# Confidential and proprietary. No part of this file may be redistributed
# without express written permission of Gimpel Software LLC.
#
# This file is provided by Gimpel Software LLC (https://www.gimpel.com) for
# use exclusively with PC-lint Plus. Redistribution to and use by licensed
# users is permitted. Any such redistribution must preserve this notice and,
# if the redistributed file has been modified, provide notice that the file
# has been modified from the original.
##############################################################################

from __future__ import print_function
import regex
import yaml
import os
import sys
import subprocess
import argparse
import string
import stat
import tempfile
import ntpath
from datetime import datetime

__version__ = "1.3.0"

def emit_note(text):
    sys.stderr.write("Note: " + text)

def emit_warning(text):
    sys.stderr.write("Warning: " + text)

def emit_error(text):
    sys.stderr.write("Error: " + text)
    sys.exit(1)

def makeHeaderGuardName(fname):
    prefix = "GS_PCLP_"
    fname = ntpath.basename(fname)
    fname = regex.sub("[.-]", "_", fname)
    fname = regex.sub("[^[:alnum:]_]", "", fname)
    return prefix + fname.upper()

def processConfig(config_file):
    try:
        config = yaml.load(open(config_file), Loader=yaml.Loader);
        return config
    except yaml.YAMLError as exc:
        emit_error("unable to parse configuration file '" + config_file + "': " + str(exc) + "\n")
    except IOError as exc:
        if not os.path.isabs(config_file):
            # If we didn't find a config file specified with
            # a relative path in the working directory, try
            # again in this script's own directory.
            script_dir = os.path.dirname(os.path.abspath(__file__))
            abs_config = os.path.join(script_dir, config_file)
            try:
                config = yaml.load(open(abs_config), Loader=yaml.Loader);
                return config
            except yaml.YAMLError as exc_abs:
                emit_error("unable to parse configuration file '" + abs_config + "': " + str(exc_abs) + "\n")
            except IOError as exc_abs:
                emit_error("unable to open configuration file '" + config_file + "' in the working directory, '" +
                 os.getcwd() + "', nor in the script directory, '" + script_dir + "': " + str(exc) + ", " + str(exc_abs) + "\n")
        emit_error("unable to open configuration file '" + config_file + "': " + str(exc) + "\n")

def runCommand(command, prog_input=None):
    # Run a command and return the collected stdout, stderr, and return value
    # Command should be list containing the executable and any arguments
    # prog_input, if provided, is sent to the stdin stream of the program.
    try:
        if prog_input is None:
            child_process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            std_output, std_error = child_process.communicate()
            exit_code = child_process.returncode
            return std_output.decode('utf-8'), std_error.decode('utf-8'), exit_code
        else:
            child_process = subprocess.Popen(command, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
            std_output, std_error = child_process.communicate(input=prog_input.encode())
            exit_code = child_process.returncode
            return std_output.decode('utf-8'), std_error.decode('utf-8'), exit_code
    except OSError as exc:
        emit_error("unable to execute command '" + " ".join(command) + "': " + str(exc) + "\n")

def listSupportedCompilers(config):
    # Dump the supported compilers and their descriptions to stdout.
    print("{:20}{}".format("Compiler Family", "Description"))
    print("{:20}{}".format("---------------", "-----------"))
    if 'compilers' in config:
        for compiler in sorted(config['compilers']):
            if 'description' in config['compilers'][compiler]:
                #print(compiler, "-", config['compilers'][compiler]['description'])
                print("{:20}{}".format(compiler, config['compilers'][compiler]['description']))


def getCompilerVersion(config, compiler, exe):
    # Run the compiler and extract the version information from it
    version_instructions = config.get('compilers', {}).get(compiler, {}).get('version', {})
    if version_instructions is None:
        emit_warning("don't know how to extract compiler version for this compiler\n")
        return None

    version = None
    if 'command' in version_instructions:
        if exe is not None:
            # We need to launch the compiler to extract the version number
            command = version_instructions['command']
            command.insert(0, exe)
            std_out, std_err, ret_val = runCommand(command)
            if 'match_expr' in version_instructions:
                result = regex.search(version_instructions['match_expr'], std_out if version_instructions.get('channel', '') == 'stdout' else std_err, regex.MULTILINE)
                if result is None:
                    emit_warning("unable to extract compiler version\n")
                    return None
                return result.group('version')
        else:
            emit_warning("need to specify compiler location with --compiler-bin to extract version information\n")
            return None

    return None

def createTemporaryFileWithContents(contents):
    t = tempfile.NamedTemporaryFile(mode='w', delete=False)
    t.file.write(contents)
    return t.name

def generateSizeOptions(config, args):
    # Find the instructions for generating PCLP size options for compiler
    # from the configuration database and execute the instructions returning
    # the result.
    compiler = args.compiler
    base_options = []
    if args.compiler_options:
        base_options = args.compiler_options.split()
    exe = args.compiler_bin
    size_instructions = config.get('compilers', {}).get(compiler, {}).get('size_options')
    if size_instructions is None:
        emit_note("size options for this compiler cannot be determined automatically, please set size options manually in the generated .lnt file\n")
        return None

    # The presence of 'command' means we should try to extract the information by
    # invoking the compiler.
    if 'command' in size_instructions and exe is not None:
        # We need to launch the compiler to extract the size options
        # The options we pass to the compiler (if any) are stored in 'command'
        # The input sent to the compiler's stdin stream, if any, is in 'input'
        command = [exe]

        if base_options:
            command = command + base_options
        if size_instructions['command']:
            command = command + size_instructions['command']

        tempfilename = None
        if 'tempfile' in size_instructions:
            tempfilename = createTemporaryFileWithContents(size_instructions['tempfile'])
            command.append(tempfilename)

        compiler_input = size_instructions.get('input')
        std_out, std_err, ret_val = runCommand(command, compiler_input)

        if tempfilename:
            os.remove(tempfilename)

        # The 'channel' determines where in the output we should look for the size options
        result_channel = size_instructions.get('channel')
        result_text = std_err if result_channel == 'stderr' else std_out

        # 'match_expr' is the pattern to use to extract the size options from
        # the output channel.  If this doesn't exist, we'll just return it all.
        match_expr = size_instructions.get('match_expr')
        if match_expr is None:
            return result_text

        match_result = regex.search(match_expr, result_text, regex.MULTILINE | regex.DOTALL)
        if match_result is not None:
            if 'size_options' in match_result.groupdict():
                # The size options should be in a named capture group called 'size_options'
                matched_portion = match_result.group('size_options')
                if matched_portion is not None:
                    return matched_portion

            if {'size_name', 'size_value'} <= set(match_result.groupdict()):
                matched_names = match_result.captures('size_name')
                matched_values = match_result.captures('size_value')

                if len(matched_names) == len(matched_values):
                    size_option_list = zip(matched_names, matched_values)
                    size_options = []
                    for size_name, size_value in size_option_list:
                        size_options.append('-' + size_name + size_value)
                    return " ".join(size_options)

    if 'fallback_values' in size_instructions:
        return size_instructions.get('fallback_values')

    if 'command' in size_instructions and exe is None:
        emit_warning("size options could not be extracted because the compiler binary was not provided, please set size options manually in the generated .lnt file or rerun with the --compiler-bin option\n")

    return None


def generateIncludeOptions(config, args):
    # Find the built-in compiler include paths, typically used to find standard
    # and system headers.  These can reside in 'cpp_include_paths', 'c_include_paths',
    # and 'include_paths'.  Well process all we find, dedupe them, and put them in
    # the order listed above.
    compiler = args.compiler
    base_options = args.compiler_options.split()
    exe = args.compiler_bin

    include_paths = config.get('compilers', {}).get(compiler, {}).get('include_paths')
    c_include_paths = config.get('compilers', {}).get(compiler, {}).get('c_include_paths')
    cpp_include_paths = config.get('compilers', {}).get(compiler, {}).get('cpp_include_paths')

    all_include_paths = []
    if cpp_include_paths is not None:
        all_include_paths.append(cpp_include_paths)
    if c_include_paths is not None:
        all_include_paths.append(c_include_paths)
    if include_paths is not None:
        all_include_paths.append(include_paths)
    if len(all_include_paths) == 0:
        return None

    found_paths = []

    for include_path in all_include_paths:
        result_text = ''
        if 'command' in include_path:
            if not exe:
                emit_warning("unable to extract include path information from compiler, use --compiler-bin to specify compiler locations\n")
                return None
            # We'll be invoking the compiler to paths
            command = [exe]
            if base_options:
                command.extend(base_options)
            command.extend(include_path['command'])
            compiler_input = include_path.get('input')
            std_out, std_err, ret_val = runCommand(command, compiler_input)

            # The 'channel' determines where in the output we should look for the include paths
            result_channel = include_path.get('channel')
            result_text = std_err if result_channel == 'stderr' else std_out
        elif 'env_var' in include_path:
            result_text = os.environ.get(include_path['env_var'], '')

        # 'match_expr' is the pattern to use to extract the include paths from
        # the output channel.  If this doesn't exist, we'll just continue.
        match_expr = include_path.get('match_expr')
        if match_expr is None:
            continue
        match_result = regex.search(match_expr, result_text, regex.MULTILINE | regex.DOTALL)
        if match_result is not None:
            # The paths should be in one or more named capture groups called 'include_dir'
            matched_portions = match_result.captures('include_dir')
            if matched_portions is not None:
                for matched_portion in matched_portions:
                    matched_portion = matched_portion.strip()
                    matched_portion = regex.sub(r'\s+', " ", matched_portion)
                    if matched_portion not in found_paths:
                        found_paths.append(matched_portion)

    return found_paths


def defaultIgnoredMacroNames():
    ignored_names = ( '_Pragma', '__BASE_FILE__', '__COUNTER__', '__DATE__',
        '__FILE__', '__INCLUDE_LEVEL__', '__LINE__', '__TIMESTAMP__',
        '__TIME__', '__VA_ARGS__', '__cplusplus', '__has_attribute', '__has_builtin',
        '__has_extension', '__has_feature', '__has_include',
        '__has_include_next', '__has_warning', '__is_identifier', 'and',
        'and_eq', 'bitand', 'bitor', 'compl', 'define', 'defined', 'not',
        'not_eq', 'or', 'or_eq', 'xor', 'xor_eq')
    return ignored_names


def extractPotentialMacros(filename):
    "Extract and return a uniq list of potential macro names from file."
    ignored_names = defaultIgnoredMacroNames()
    names = set()
    try:
        with open(filename, "rb") as f:
            s = ""
            for c in f.read():
                if (c in string.letters + '_') or (s and c in string.digits):
                    s += c
                    continue
                elif s:
                    if s not in ignored_names:
                        names.add(s)
                    s = ""
    except IOError as exc:
        emit_error("unable to scavenge macros from file '" + filename + "': " + str(exc) + "\n")
    return names


def createScavengeData(names):
    # Produce scavenger data for preprocessing
    scavenge_data = ''
    for name in names:
        scavenge_data += ("#ifdef %s\n-d%s{%s}\n#endif\n" % (name, name, name))
    return scavenge_data


def extractScavengeResults(output):
    # Remove anything except -d options from preprocessed scavenger output.
    good_lines = []
    good_pattern = regex.compile("^\s*-d")

    for line in output.split('\n'):
        if good_pattern.match(line):
            # Convert to #define
            line = regex.sub("^\s*-d(.*?)\{(.*)\}$", "#define \\1 \\2", line)
            good_lines.append(line)

    return good_lines


def generateMacroDefinitions(instructions, args, base_options):
    # There are several ways macro definitions may be generated:
    #   'command' - via compiler invocation and pattern matching
    #   'definitions' - explicitly specified
    #   'scavenge' - using the macro scavenger method
    #
    # For 'command', the macro definitions are expected to appear in the
    # output of the invocation, either in 'stdout' or 'stderr' as indicated
    # by 'channel'.  'match_expr' will be used to match the definitions, if
    # present, otherwise the output will be used as-is.
    #
    # For 'definitions', a list of lists is expected where the first item
    # in each list is the definition name, including parameter list for
    # function-like macros, and the second item is the definition, or null
    # for no definition.  E.g.:
    #
    #   [['A', 1], ['B', ''], ['C', null], ['f(a,b)', '(a + b)']]
    #
    # results in the definitions:
    #
    #   #define A 1
    #   #define B
    #   #define C
    #   #define f(a,b) (a + b)
    #
    # For 'scavenge', 'scavenge_files' and/or 'scavenge_dirs', provided
    # on the command line, are used to search for possible macros.
    # A scavenger file is built and passed to the compiler using the
    # preprocessor to expand defined macros and the output is exhumed
    # to dig out and build a macro list.
    #
    # In all cases, a list of macro definitions, one per line, is returned.
    exe = args.compiler_bin
    scavenging = args.scavenge_files or args.scavenge_dirs

    if scavenging and 'scavenge_command' in instructions:
        scavenge_files = set()
        if args.scavenge_files:
            scavenge_files.update(args.scavenge_files)
        if args.scavenge_dirs:
            for scavenge_dir in args.scavenge_dirs:
                for folder, subs, files in os.walk(scavenge_dir):
                    for filename in files:
                        if args.scavenge_pattern:
                            if not regex.match(args.scavenge_pattern, filename):
                                continue
                        full_path = os.path.join(folder, filename)
                        stat_info = os.stat(full_path)
                        if stat.S_ISREG(stat_info.st_mode):
                            scavenge_files.add(full_path)

        potential_macros = set()
        for scavenge_file in scavenge_files:
            potential_macros.update(extractPotentialMacros(scavenge_file))

        scavenge_data = createScavengeData(potential_macros)

        command = [exe] + instructions['scavenge_command']
        compiler_input = scavenge_data
        std_out, std_err, ret_val = runCommand(command, compiler_input)

        result_channel = instructions.get('channel')
        result_text = std_err if result_channel == 'stderr' else std_out

        return "\n".join(extractScavengeResults(result_text)) + "\n"

    if 'command' in instructions:
        if exe is None:
            emit_warning("unable to extract macro definitions from compiler, use --compiler-bin to specify compiler locations\n")
        else:
            command = [exe]
            if base_options:
                command.extend(base_options)
            command.extend(instructions['command'])

            tempfilename = None
            if 'tempfile' in instructions:
                tempfilename = createTemporaryFileWithContents(instructions['tempfile'])
                command.append(tempfilename)

            compiler_input = instructions.get('input')
            std_out, std_err, ret_val = runCommand(command, compiler_input)

            if tempfilename:
                os.remove(tempfilename)

            result_channel = instructions.get('channel')
            result_text = std_err if result_channel == 'stderr' else std_out

            # 'match_expr' is the pattern to use to extract the macros  from
            # the output channel.
            match_expr = instructions.get('match_expr')
            if match_expr is None:
                return result_text

            match_result = regex.search(match_expr, result_text, regex.MULTILINE|regex.DOTALL)
            if match_result:
                ignored_names = defaultIgnoredMacroNames()
                macros = match_result.group('macros')
                macros_to_keep = []
                for md in macros.splitlines():
                    m = regex.match(r'\s*#define\s+(\w+)', md)
                    if not m or m.group(1) in ignored_names:
                        continue
                    macros_to_keep.append(md)
                return "\n".join(macros_to_keep) + "\n"

    if 'definitions' in instructions:
        macro_defs_str = ''
        for definition in instructions['definitions']:
            macro_name, macro_def = definition
            macro_defs_str += '#define ' + macro_name + ' ' + macro_def + "\n"
        return macro_defs_str

    return None

def generateDecls(config, args):
    # Find the built-in compiler decls.  We can have 'c_decls', 'cpp_decls',
    # and 'decls'.
    compiler = args.compiler
    base_options = args.compiler_options.split()
    exe = args.compiler_bin

    decls = config.get('compilers', {}).get(compiler, {}).get('decls', {}).get('definitions')
    c_decls = config.get('compilers', {}).get(compiler, {}).get('c_decls', {}).get('definitions')
    cpp_decls = config.get('compilers', {}).get(compiler, {}).get('cpp_decls', {}).get('definitions')

    all_decls = ''
    if decls:
        all_decls += "\n".join(decls)
    if c_decls:
        all_decls += '#ifndef __cplusplus\n' + "\n".join(c_decls) + '\n#endif\n'
    if cpp_decls:
        all_decls += '#ifdef __cplusplus\n' + "\n".join(cpp_decls) + '\n#endif\n'

    return all_decls


def generateMacros(config, args):
    # Find the built-in compiler macros.  We can have 'c_macros', 'cpp_macros',
    # and 'macros'.
    compiler = args.compiler
    base_options = args.compiler_options.split()
    exe = args.compiler_bin

    macros = config.get('compilers', {}).get(compiler, {}).get('macros')
    c_macros = config.get('compilers', {}).get(compiler, {}).get('c_macros')
    cpp_macros = config.get('compilers', {}).get(compiler, {}).get('cpp_macros')

    all_macro_definitions = ''
    if macros:
        generic_macro_defs = generateMacroDefinitions(macros, args, args.compiler_options.split())
        if generic_macro_defs:
            all_macro_definitions += generic_macro_defs
    if c_macros:
        c_macro_defs = generateMacroDefinitions(c_macros, args, args.compiler_options.split() + args.compiler_c_options.split())
        if c_macro_defs:
            all_macro_definitions += '#ifndef __cplusplus\n' + c_macro_defs + '#endif\n'
    if cpp_macros:
        cpp_macro_defs = generateMacroDefinitions(cpp_macros, args, args.compiler_options.split() + args.compiler_cpp_options.split())
        if cpp_macro_defs:
            all_macro_definitions += '#ifdef __cplusplus\n' + cpp_macro_defs + '#endif\n'

    return all_macro_definitions


def generateBaseOptions(config, args):
    compiler = args.compiler
    base_config = config.get('compilers', {}).get(compiler, {}).get('base_config', {})
    base_options = ""
    for key in sorted(base_config.keys()):
        base_options += "//     " + key.title() + "\n"
        for option, annotation in base_config[key]:
            base_options += option
            if annotation:
                base_options += "  // " + annotation
            base_options += "\n"
        base_options += "\n"

    return base_options


def generateCompilerConfig(config, args):
    if not 'compilers' in config:
        emit_error("compiler database doesn't contain any compiler data\n")

    compiler = args.compiler
    compile_commands = args.compiler_options.split()
    exe = args.compiler_bin

    if compiler is None:
        emit_error("no --compiler specified\n")
    if not compiler in config['compilers']:
        emit_error("'" + compiler + "' is not recognized for automatic configuration; use the option '--list-compilers' to list compilers in the database (refer to the 'System Configuration' section in chapter 2 of the manual to configure an unknown compiler manually)\n")

    compiler_entry = config['compilers'][compiler]
    compiler_version = getCompilerVersion(config, compiler, exe)

    intro_string = "/* Compiler configuration for %s %s.\n   Using the options: %s \n   Generated on %s with pclp_config version %s.\n */\n\n" % (compiler, compiler_version, " ".join(compile_commands), datetime.now().strftime('%Y-%m-%d %H:%M:%S'), __version__)

    # Base configuration
    base_string = ""
    base_options = generateBaseOptions(config, args)
    if base_options is None:
        base_string = "// Unable to generate base compiler options.\n\n"
    else:
        base_string = "// Base Options\n" + base_options + "\n\n"

    # Size Options
    size_string = ""
    size_options = generateSizeOptions(config, args)
    if size_options is None:
        emit_warning("unable to determine size options, these will need to be set manually in the generated .lnt file\n")
        size_string = "// Unable to determine size options. \n\n"
    else:
        size_string = "// Size Options\n" + size_options + "\n\n"

    # Built-in include directories
    includes_string = ""
    include_directories = generateIncludeOptions(config, args)
    if include_directories is None:
        emit_warning("unable to determine built-in include directories, these will need to be set manually in the generated .lnt file\n")
        includes_string = "// Failed to extract include paths. \n"
    else:
        includes_string = "// Include Options\n"
        for include_dir in include_directories:
            includes_string += "--i\"" + include_dir.strip() + "\"\n"
    includes_string += "\n"

    # Built-in macros
    builtin_macros = generateMacros(config, args)
    builtin_decls = generateDecls(config, args)

    # Custom compile commands
    custom_options_string = "// Transformed compiler options\n"
    while compile_commands:
        transformations, options_consumed = handleCompilerOption(config, compiler, compile_commands)
        if transformations:
            custom_options_string += " ".join(transformations) + " // From compiler option(s): " + " ".join(compile_commands[:options_consumed]) + "\n"
        compile_commands = compile_commands[options_consumed:]
    custom_options_string += "\n"

    if args.config_output_lnt_file:
        with open(args.config_output_lnt_file, 'w') as f:
            f.write(intro_string)
            f.write(base_string)
            f.write(custom_options_string)
            f.write(size_string)
            f.write(includes_string)
            if args.config_output_header_file:
                header_guard_macro = makeHeaderGuardName(args.config_output_header_file)
                with open(args.config_output_header_file, 'w') as h:
                    h.write('#ifndef ' + header_guard_macro + "\n")
                    h.write('#define ' + header_guard_macro + "\n")
                    h.write(builtin_macros)
                    h.write(builtin_decls)
                    h.write('\n#endif /* ' + header_guard_macro + " */\n")
                header_path = args.config_output_header_file
                if args.header_option_use_enclosing_directory:
                    header_path = "%ENCLOSING_DIRECTORY%/" + header_path
                header_string = '+libh(' + header_path + ')\n' + '-header(' + header_path + ')\n'
                f.write(header_string)
            else:
                emit_warning("no --config-output-header-file specified\n")
    else:
        emit_warning("no --config-output-lnt-file specified\n")

def generateProjectConfig(config, args):
    compiler = args.compiler
    source_pattern = args.source_pattern
    imposter_file = args.imposter_file

    if not imposter_file:
        emit_error("An imposter input file must be specified when using --generate-project-config, " +
            "use --imposter-file to specify the imposter file\n")

    imposter_contents = ''
    with open(imposter_file) as f:
        imposter_contents = '[' + ",".join(f.readlines()) + ']'

    compile_commands = yaml.load(imposter_contents, Loader=yaml.Loader)

    output_file = open(args.config_output_lnt_file, 'w')

    for compile_command in compile_commands:
        source_files = []
        options_str = ''
        while compile_command:
            sf_match = regex.match(source_pattern, compile_command[0])
            if sf_match:
                source_files.append(compile_command[0])
                compile_command.pop(0)
                continue
            transformations, options_consumed = handleCompilerOption(config, compiler, compile_command)
            if transformations:
                options_str += "\n".join(transformations) + "\n"
            compile_command = compile_command[options_consumed:]
        if source_files:
            output_file.write("-env_push\n")
            output_file.write(options_str)
            for source_file in source_files:
                output_file.write('"' + source_file + '"\n')
            output_file.write("-env_pop\n\n")

    output_file.close()


def handleCompilerOption(config, compiler, option_list):
    # Given a configuration and a compiler, attempt to decode the provided
    # compiler option, mapping it to the corresponding PC-lint option if
    # an appropriate transformation exists.  Returns the transformation
    # (or None) and the number of options consumed.
    option_map = config.get('compilers', {}).get(compiler, {}).get('options')
    if option_map is None or not option_list:
        return (None, 1)
    option_str = option_list[0]

    # Find the longest matching option
    best_match = None
    best_match_size = 0
    for candidate_option in option_map:
        if option_str.startswith(candidate_option):
            if len(candidate_option) > best_match_size:
                best_match = candidate_option
                best_match_size = len(candidate_option)

    if best_match is None:
        # We didn't recognize the option
        return (None, 1)

    found_option_map = option_map[best_match]
    if found_option_map is None:
        return (None, 1)

    if 'transform' in found_option_map:
        return ([found_option_map['transform']], 1)

    if 'transforms' in found_option_map:
        options_consumed = 1
        while True:
            replacements = []
            for transform_pair in found_option_map['transforms']:
                match_pattern, repl_pattern = transform_pair
                replaced_option_str = option_str
                if 'pre_transforms_replacements' in found_option_map:
                    for replacement_pair in found_option_map['pre_transforms_replacements']:
                        replaced_option_str = replaced_option_str.replace(replacement_pair[0], replacement_pair[1])
                match_result = regex.match(match_pattern, replaced_option_str)
                if match_result:
                    replacement = regex.sub(match_pattern, repl_pattern, replaced_option_str)
                    replacements.append(replacement)
            if replacements:
                return (replacements, options_consumed)
            # Didn't match any of the transformation patterns.  Add the next
            # option to the option string to see if we can match with an arg.
            if len(option_list) > options_consumed:
                option_str += " " + option_list[options_consumed]
                options_consumed += 1
            else:
                return (None, 1)

    # Found an option but no transformations exist
    return (None, 1)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--list-compilers',
        help='list the supported compilers',
        default=False,
        action='store_true')

    parser.add_argument('--generate-compiler-config',
        help='generate a customized compiler configuration',
        default=False,
        action='store_true')

    parser.add_argument('--generate-project-config',
        help='generate a customized project configuration',
        default=False,
        action='store_true')

    parser.add_argument('--compiler',
        help='the compiler that will be used to generate configurations',
        type=str)

    parser.add_argument('--compiler-bin',
        help='the location of the compiler executable',
        type=str)

    parser.add_argument('--compiler-options',
        help='base compiler options',
        default='')

    parser.add_argument('--compiler-c-options',
        help='base C-specific compiler options',
        default='')

    parser.add_argument('--compiler-cpp-options',
        help='base C++-specific compiler options',
        default='')

    parser.add_argument('--ignore-options',
        help='compiler options that should be ignored (not transformed)',
        nargs='+')

    parser.add_argument('--compiler-database',
        help='the name of the compiler database file',
        default='compilers.yaml')

    parser.add_argument('--repl',
        help="enter compiler options and see the transformations that would be made",
        default=False,
        action='store_true')

    parser.add_argument('--compiler-version',
        help="show the version of the compiler being configured",
        default=False,
        action='store_true')

    parser.add_argument('--source-pattern',
        help="the pattern used to match project source files in compiler invocations",
        default='.*\.(c|cpp)$',
        type=str)

    parser.add_argument('--imposter-file',
        help="the file containing compiler invocations logged by the imposter program",
        type=str)

    parser.add_argument('--config-output-lnt-file',
        help="the file to write the configuration to",
        type=str)

    parser.add_argument('--config-output-header-file',
        help="the file to write supplemental configuration data to (macro definitions, etc)",
        type=str)

    parser.add_argument('--scavenge-files',
        help="the list of files to attempt to extract macro information from",
        action='append')

    parser.add_argument('--scavenge-dirs',
        help="the list of directories to recursively process files from to extract macro information from",
        action='append')

    parser.add_argument('--scavenge-pattern',
        help="the regular expression pattern used to match filenames, excluding path, for macro extraction",
        type=str)

    parser.add_argument('--header-option-use-enclosing-directory',
        help="use the built-in %%ENCLOSING_DIRECTORY%% environment variable to provide an 'absolute' path for the compiler configuration -header option",
        default=False,
        action='store_true')

    args = parser.parse_args()

    handled_task = False

    if args.compiler_version:
        config = processConfig(args.compiler_database)
        print(getCompilerVersion(config, args.compiler, args.compiler_bin))
        handled_task = True

    if args.list_compilers:
        config = processConfig(args.compiler_database)
        listSupportedCompilers(config)
        handled_task = True

    if args.generate_compiler_config:
        config = processConfig(args.compiler_database)
        generateCompilerConfig(config, args)
        handled_task = True

    if args.generate_project_config:
        config = processConfig(args.compiler_database)
        generateProjectConfig(config, args)
        handled_task = True

    if args.repl:
        config = processConfig(args.compiler_database)
        handled_task = True
        while True:
            line = sys.stdin.readline()
            if not line:
                break
            print(handleCompilerOption(config, args.compiler, line.strip().split()))

    if not handled_task:
        if args.imposter_file:
            emit_warning("No work done as no task was requested, did you forget the --generate-project-config option?\n")
        elif args.config_output_lnt_file or args.config_output_header_file:
            emit_warning("No work done as no task was requested, did you forget the --generate-compiler-config option?\n")
        else:
            emit_warning("No work done as no task was requested, use --help for usage.\n")

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        emit_error(str(e) + "\n")
