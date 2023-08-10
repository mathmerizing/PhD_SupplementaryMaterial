import os

# read in all log files from logs/
n = 0
labels = []
for file in os.listdir('logs/'):
    # get r and s from file name
    r = int(file.split('_')[1].strip('r'))
    s = int(file.split('_')[2].split('.')[0].strip('s'))
    print(f"r = {r}, s = {s}")
    if not(
        r == 0 and s == 1
        or r == 0 and s == 2
        or r == 1 and s == 1
        or r == 1 and s == 2
        or r == 2 and s == 1
        or r == 2 and s == 2
        or r == 2 and s == 3
        or r == 3 and s == 2
        or r == 3 and s == 3
        or r == 3 and s == 4
    ):
        print("Skipping")
        continue
    
    if file.endswith(".txt"):
        with open('logs/'+file) as f:
            # read in all lines from log file
            lines = f.readlines()

            # get the last lines of the log file, which start with cycle
            last_lines = []
            is_conv_table = False
            for line in lines:
                if line.startswith('cycle'):
                    is_conv_table = True
                if line.startswith('n'):
                    is_conv_table = False
                if is_conv_table:
                    last_lines.append(line)

            # for l in last_lines:
            #     print(l)

            dofs = []
            error = []
            for l in last_lines[1:]:
                if l == "\n":
                    continue
                #print(f"l = '{l}'")
                dofs.append(int(l.split()[2]))
                error.append(float(l.split()[5]))

            print(f"\\addplot[color=brewer{n%6+1},mark={'otimes' if n < 6 else 'diamond'},style=ultra thick] coordinates "  + "{", end="")
            for i in range(len(dofs)):
                print(f"({dofs[i]},{error[i]})", end="")
            print("};")

            labels.append(f"$\\cG({s})\\dG({r})$")
            n += 1
            
            #quit()
print("\n\n" + "\\legend{" + ",".join(labels) + "}")