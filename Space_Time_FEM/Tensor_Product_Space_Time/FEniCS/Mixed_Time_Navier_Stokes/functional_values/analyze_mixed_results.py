import os
import numpy as np
import matplotlib.pyplot as plt

# for each directory in results_mixed_dG/2D-3 compute the error in drag and lift
lobatto = {"drag": {}, "lift": {}}
legendre = {"drag": {}, "lift": {}}

# load drag_featflow.dat and lift_featflow.dat files
time_reference, drag_reference = [], []
with open('featflow/drag_featflow.dat', 'r') as f:
    lines = f.readlines()
    for line in lines[1:]:
        line = line.split(',')
        time_reference.append(float(line[0]))
        drag_reference.append(float(line[1]))

# integrate drag_reference in time
drag_qoi_reference = np.trapz(time_reference, drag_reference) / 8.

lift_reference = []
with open('featflow/lift_featflow.dat', 'r') as f:
    lines = f.readlines()
    for line in lines[1:]:
        line = line.split(',')
        lift_reference.append(float(line[1]))

# integrate lift_reference in time
lift_qoi_reference = np.trapz(time_reference, lift_reference) / 8.

# plot reference drag and lift values
# plt.figure()
# plt.plot(time_reference, drag_reference)
# plt.xlabel('Time [s]')
# plt.ylabel('Drag [N]')
# plt.title('Drag over time for reference')
# plt.show()

# plt.figure()
# plt.plot(time_reference, lift_reference)
# plt.xlabel('Time [s]')
# plt.ylabel('Lift [N]')
# plt.title('Lift over time for reference')
# plt.show()


for folder in os.listdir('results_mixed_dG/2D-3'):
    #print(folder)

    # load drag_values.npy, lift_values.npy, and times_draglift.npy
    drag_values = np.load('results_mixed_dG/2D-3/' + folder + '/drag_values.npy')
    lift_values = np.load('results_mixed_dG/2D-3/' + folder + '/lift_values.npy')
    times_draglift = np.load('results_mixed_dG/2D-3/' + folder + '/times_draglift.npy')

    # concatenate drag_values and lift_values and times_draglift into one array
    values = np.vstack((times_draglift, drag_values, lift_values)).transpose()
    # print(values.shape)
    
    # remove each row of values which has a nan value
    values = values[~np.isposinf(values).any(axis=1)]
    # print(values.shape)

    # split values into times_draglift, drag_values, and lift_values
    times_draglift = values[:,0]
    drag_values = values[:,1]
    lift_values = values[:,2]

    # integrate drag_values in time
    drag_qoi = np.trapz(times_draglift, drag_values) / 8.

    # integrate lift_values in time
    lift_qoi = np.trapz(times_draglift, lift_values) / 8.

    # compute error in drag and lift
    drag_error = abs(drag_qoi - drag_qoi_reference)
    lift_error = abs(lift_qoi - lift_qoi_reference)

    # store error in drag and lift in lobatto or legendre
    if "Lobatto" in folder:
        lobatto["drag"][folder] = drag_error
        lobatto["lift"][folder] = lift_error
    else:
        legendre["drag"][folder] = drag_error
        legendre["lift"][folder] = lift_error

    if "Lobatto" in folder:
        # create drag.dat and lift.dat files
        r_v = int(folder.split('_')[2])
        r_p = int(folder.split('_')[4])
        
        # create drag.dat and lift.dat files

        # lift
        with open('featflow/lift.dat', 'w') as f:
            f.write('x,y\n')
            for i in range(len(times_draglift)):
                f.write(str(times_draglift[i]) + ',' + str(lift_values[i]) + '\n')
        
        # drag
        with open('featflow/drag.dat', 'w') as f:
            f.write('x,y\n')
            for i in range(len(times_draglift)):
                f.write(str(times_draglift[i]) + ',' + str(drag_values[i]) + '\n')

        # change into featflow directory
        os.chdir('featflow')

        # create a backup of figure.tex
        os.system('cp figure.tex figure_backup.tex')

        # replace the METHODNAME in figure.tex
        METHODNAME = "$" + ' r_v = ' + str(r_v) + '| r_p = ' + str(r_p) + "$"
        os.system('sed -i "s/METHODNAME/' + METHODNAME + '/g" figure.tex')

        # compile figure.tex
        os.system('bash make_figures.sh')

        # assert that lift.pdf and drag.pdf exist
        assert os.path.isfile('lift.pdf'), "lift.pdf does not exist"
        assert os.path.isfile('drag.pdf'), "drag.pdf does not exist"

        # move lift.pdf and drag.pdf to results_mixed_dG/2D-3/folder
        os.system('mv lift.pdf ../results_mixed_dG/2D-3/' + folder + '/')
        os.system('mv drag.pdf ../results_mixed_dG/2D-3/' + folder + '/')

        # remove figure.tex
        os.system('rm figure.tex')

        # restore figure_backup.tex
        os.system('mv figure_backup.tex figure.tex')

        # change back to main directory
        os.chdir('..')

print(lobatto)

errors = {quad: {qoi: -1. * np.ones(((3,3))) for qoi in ["drag", "lift"]} for quad in ["Lobatto", "Legendre"]}
for name, table in zip(["Lobatto", "Legendre"], [lobatto, legendre]):
    print(name)

    for qoi in ["drag", "lift"]:
        print(qoi)
        
        for key in table[qoi]:
            # get r_v and r_p from key
            r_v = int(key.split('_')[2])
            r_p = int(key.split('_')[4])
            errors[name][qoi][r_v, r_p] = table[qoi][key]

        print(errors[name][qoi])

# print(errors)

# create a latex table with the errors in drag and lift
# print("Lobatto")
# print("drag")
# for key in lobatto["drag"]:
#     print(key, "&", lobatto["drag"][key], "\\\\")
# print("lift")
# for key in lobatto["lift"]: