def parse_out(filenames):
    n = len(filenames)
    files = []
    for i in range(n):
        file = open(filenames[i],'r')
        files.append(file.readlines())

    for i in range(n):
        setting_name = (filenames[i].split('.'))[-2]
        with open(setting_name + '_affected.txt', 'a') as the_file:
            the_file.write(setting_name + '\n' + '----------------------\n' )
        for j in range(len(files[i])):
            line = files[i][j]
            line = line.strip()

            if line.startswith('read problem'):
                line = line.split('/')
                name = line[-1]
                name = name[:-9]
                next_instance = True

            elif line.startswith('Quadratic Nlhdlr :'):
                line = line.split(' ')
                line = filter(None,line)

                if (float(line[3]) > 0) and (next_instance == True):
                    print(name)
                    with open(setting_name + '_affected.txt', 'a') as the_file:
                        the_file.write(name + '\n')
                    next_instance = False
parse_out(['check/results/check.MINLP_quadratics.scip.opt.strengthening_nocutlimit.out'])
