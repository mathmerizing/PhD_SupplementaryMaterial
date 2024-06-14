import os
import numpy as np
import matplotlib.pyplot as plt

# for each directory in results_equal_dG create drag and lift plots
for folder in os.listdir('results_equal_dG/2D-3'):
    METHODNAME = ""
    if "fine" in folder:
        METHODNAME = "fine mesh"
    elif "pygmsh" in folder:
        METHODNAME = "pygmsh mesh"
    else:
        METHODNAME = "coarse mesh"
    
    # load drag_values.npy, lift_values.npy, and times_draglift.npy
    drag_values = np.load('results_equal_dG/2D-3/' + folder + '/drag_values.npy')
    lift_values = np.load('results_equal_dG/2D-3/' + folder + '/lift_values.npy')
    times_draglift = np.load('results_equal_dG/2D-3/' + folder + '/times_draglift.npy')
    
    # concatenate drag_values and lift_values and times_draglift into one array
    values = np.vstack((times_draglift, drag_values, lift_values)).transpose()
    # print(values.shape)
    
    # remove each row of values which has a nan value
    values = values[~np.isposinf(values).any(axis=1)]
    # print(values.shape)

    # remove the number of rows because TeX can't handle too many data points
    values = values[::3, :]
    
    # split values into times_draglift, drag_values, and lift_values
    times_draglift = values[:,0]
    drag_values = values[:,1]
    lift_values = values[:,2]
    
    # plot drag and lift values

    # drag

    # plt.figure()
    # plt.plot(times_draglift, drag_values)
    # plt.xlabel('Time [s]')
    # plt.ylabel('Drag [N]')
    # plt.title('Drag over time for ' + METHODNAME)
    # plt.show()

    # # lift
    # plt.figure()
    # plt.plot(times_draglift, lift_values)
    # plt.xlabel('Time [s]')
    # plt.ylabel('Lift [N]')
    # plt.title('Lift over time for ' + METHODNAME)
    # plt.show()

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
    
    # change into feaflow directory
    os.chdir('featflow')

    # create a backup of figure.tex
    os.system('cp figure.tex figure_backup.tex')

    # replace the METHODNAME in figure.tex
    os.system('sed -i "s/METHODNAME/' + METHODNAME + '/g" figure.tex')

    # compile figure.tex
    os.system('bash make_figures.sh')

    # move lift.pdf and drag.pdf to results_equal_dG/2D-3/folder
    os.system('mv lift.pdf ../results_equal_dG/2D-3/' + folder + '/')
    os.system('mv drag.pdf ../results_equal_dG/2D-3/' + folder + '/')

    # remove figure.tex
    os.system('rm figure.tex')

    # restore figure_backup.tex
    os.system('mv figure_backup.tex figure.tex')

    # change back to main directory
    os.chdir('..')