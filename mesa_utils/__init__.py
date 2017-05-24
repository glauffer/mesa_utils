import mesa_reader as ms
import matplotlib.pyplot as plt
import numpy as np

__all__ = ['hrd', 'abun_plot']


def hrd(folder, title=None, save=False, name=None):
    '''
    Function to plot the Hertzsprung-Russel diagram

    Parameter
    ---------
    folder  : str
            path to LOGS folder -> '/path/to/LOGS/'
    title   : str
            title to use on the plot -> 'title'
    save    : bool
            if True, save a .png file of the plot
    name    : str
            name of the .png file if the parameter save is True

    Returns
    -------
    show the plot or save a .png file
    '''

    path = ms.MesaLogDir(folder)
    hist = path.history
    lum = hist.data('log_L')
    teff = hist.data('log_Teff')
    plt.figure(figsize=(10, 10))
    plt.plot(teff, lum)
    plt.gca().invert_xaxis()
    plt.suptitle('HR Diagram', fontsize=18)
    plt.title(title)
    plt.xlabel(r'$\log T_{Eff}$', fontsize=16)
    plt.ylabel(r'$\log L$', fontsize=16)

    if save is True:
        plt.savefig(name + '.png')
        return 'File save with name ' + name
    else:
        plt.show()


def abun_plot(folder, mod_n=None, x_lim=12, title=None,
              save=False, name=None, isotope=None):
    '''
    Plot the abundance profile at model number

    Parameter
    ---------
    folder  : str
            path to LOGS folder -> '/path/to/LOGS/'
    mod_n   : int
            model number to calculate the profiles
    x_lim   : int
            limit to x-axis on profile plot. Default = 12
    title   : str
            title to use on the plot -> 'title'
    save    : bool
            if True, save a .png file of the plot
    name    : str
            name of the .png file if the parameter save is True
    isotope: Not implemented yet

    Returns
    -------
    show the plot or save a .png file
    '''
    # Talvez implementar um parametro para escolher quais isotopos
    # serao plotados
    # if isotope == 'CO':
    #     isos = ['h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24']
    # elif isotope == 'ONe':
    #     isos = ['h1', 'he3','he4', 'c12', 'n14', 'o16', 'ne20', 'mg24']

    path = ms.MesaLogDir(folder)
    hist = path.history
    # models = hist.data('model_number')
    mass = hist.data('star_mass')
    teff = hist.data('log_Teff')

    prof = path.profile_data(mod_n)
    # isos = ['h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24']
    # for i in isos:
    x = - prof.data('logxq')

    h1 = prof.data('h1')
    he4 = prof.data('he4')
    c12 = prof.data('c12')
    n14 = prof.data('n14')
    o16 = prof.data('o16')
    ne20 = prof.data('ne20')
    mg24 = prof.data('mg24')

    plt.figure(figsize=(14, 7))
    plt.plot(x, he4, color='darkolivegreen', lw=1.5, label=r'$He4$')
    plt.plot(x, c12, color='k', lw=1.5, label=r'$C12$')
    plt.plot(x, n14, color='r', lw=1.5, label=r'$N14$')
    plt.plot(x, o16, color='g', lw=1.5, label=r'$O16$')
    plt.plot(x, ne20, color='y', lw=1.5, label=r'$Ne20$')
    plt.plot(x, mg24, color='m', lw=1.5, label=r'$Mg24$')
    plt.plot(x, h1, color='b', lw=1.5, label=r'$H1$')

    plt.xlim(0, x_lim)
    plt.xlabel(r'$\log (1 - q)$')
    plt.ylabel('mass fraction')
    plt.suptitle(title)
    plt.title('Abundance at model ' + str(mod_n) + r' $M_i = $ ' +
              str(mass[0]) + r' $M_f = $ ' + str(mass[mod_n]) +
              r' $Teff = $ ' + str(10**teff[mod_n]))

    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

    if save is True:
        plt.savefig(name + '.png')
        return 'File save with name ' + name
    else:
        plt.show()


def get_centerMass(path):
    '''
    Extract the central mass of h1, he4, c12, n14, o16, ne20 and mg24
    from the LOGS directory

    Parameters
    ----------
    path    : string
            LOGS directory path

    Returns
    -------
    numpy array with initial mass and central values of h1, he4, c12, n14,
    o16, ne20 and mg24
    '''

    fol = ms.MesaLogDir(log_path=path)
    h = fol.history
    c_h1 = h.data('center_h1')
    c_he4 = h.data('center_he4')
    c_c12 = h.data('center_c12')
    c_n14 = h.data('center_n14')
    c_o16 = h.data('center_o16')
    c_ne20 = h.data('center_ne20')
    mi = h.data('star_mass')

    prof_numbers = fol.profile_numbers
    prof = fol.profile_data(profile_number=prof_numbers[-1])
    c_mg24 = prof.data('mg24')

    return np.array([mi[0], c_h1[-1], c_he4[-1], c_c12[-1],
                     c_n14[-1], c_o16[-1], c_ne20[-1], c_mg24[0]])
    # return np.array([mass, c_c12[-1], c_o16[-1], c_ne20[-1], c_mg24[0]])
    # return [mass, ' | ', c_c12[-1], ' | ', c_o16[-1], ' | ', c_ne20[-1],
    #         ' | ', c_mg24[-1]]


def print_c_mass(list_path):
    '''
    Print the initial mass and center values of h1, he4, c12, n14, o16, ne20
    and mg24

    Parameter
    --------
    list_path   : list of strings
            list of all LOGS directories that are going to be analized
    '''

    np.set_printoptions(formatter={'float': '{: 0.5}'.format})
    # print('Mass (Msun), center c12, center o16, center ne20, center mg24')
    print('Mass(Mo), c_h1, c_he4, c_c12, c_n14, c_o16, c_ne20, c_mg24')
    for i in list_path:
        c_mass = get_centerMass(i)
        print(c_mass)
