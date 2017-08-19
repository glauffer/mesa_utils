import mesa_reader as ms
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.ndimage
import os

__all__ = ['hrd', 'abun_plot', 'get_mat_mcenter', 'print_data']


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
    mass = hist.data('star_mass')
    plt.figure(figsize=(10, 8))
    plt.plot(teff, lum)
    plt.gca().invert_xaxis()
    plt.suptitle('HR Diagram', fontsize=18)
    plt.title('Mi = ' + str(mass[0]) + ' - ' + title)
    plt.xlabel(r'$\log T_{Eff}$', fontsize=16)
    plt.ylabel(r'$\log L$', fontsize=16)

    if save is True:
        plt.savefig(name + '.pdf')
        return 'File save with name ' + name
    else:
        plt.show()


def abun_plot(folder, mod_n=None, x_lim=12, title=None,
              save=False, name=None, isotope=None, x_axis='atm'):
    '''
    Plot the abundance profile at model number using log(1-q) on x-axis
    (better to see information on atmosfere)

    Parameter
    ---------
    folder  : str
            path to LOGS folder -> '/path/to/LOGS/'
    mod_n   : int
            model number to calculate the profiles. If None, the last model
            will be used
    x_lim   : int
            limit to x-axis on profile plot. Default = 12
    title   : str
            title to use on the plot -> 'title'
    save    : bool
            if True, save a .png file of the plot
    name    : str
            name of the .png file if the parameter save is True
    x_axis  : str
            can be 'atm' or 'nuc'. If 'atm' is chosen, x-axis will be -log(1-q)
            (better for atmosfere region). If 'nuc', x-axis will be q (better
            for nucleus region). Default is 'atm'
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

    if mod_n is None:
        profiles = path.model_numbers
        model = profiles[-1]
        m_ind = -1
    else:
        model = mod_n
        m_ind = -1  # mod_n

    prof = path.profile_data(model)
    # isos = ['h1', 'he4', 'c12', 'n14', 'o16', 'ne20', 'mg24']
    # for i in isos:
    if x_axis == 'atm':
        x = - prof.data('logxq')
        xlabel = r'$\log (1 - q)$'
    else:
        x = prof.data('q')
        x_lim = 1
        xlabel = r'$\frac{m}{M}$'

    h1 = prof.data('h1')
    he3 = prof.data('he3')
    he4 = prof.data('he4')
    c12 = prof.data('c12')
    c13 = prof.data('c13')
    n13 = prof.data('n13')
    n14 = prof.data('n14')
    n15 = prof.data('n15')
    o16 = prof.data('o16')
    o18 = prof.data('o18')
    ne20 = prof.data('ne20')
    ne22 = prof.data('ne22')
    mg24 = prof.data('mg24')

    plt.figure(figsize=(14, 7))
    plt.plot(x, he3, color='c', lw=1.5, label=r'$He3$')
    plt.plot(x, he4, color='darkolivegreen', lw=1.5, label=r'$He4$')
    plt.plot(x, c12, color='k', lw=1.5, label=r'$C12$')
    plt.plot(x, c13, color='grey', lw=1.5, label=r'$C13$')
    plt.plot(x, n13, color='darkred', lw=1.5, label=r'$N13$')
    plt.plot(x, n14, color='r', lw=1.5, label=r'$N14$')
    plt.plot(x, n15, color='salmon', lw=1.5, label=r'$N15$')
    plt.plot(x, o16, color='g', lw=1.5, label=r'$O16$')
    plt.plot(x, o18, color='limegreen', lw=1.5, label=r'$O18$')
    plt.plot(x, ne20, color='y', lw=1.5, label=r'$Ne20$')
    plt.plot(x, ne22, color='orange', lw=1.5, label=r'$Ne22$')
    plt.plot(x, mg24, color='m', lw=1.5, label=r'$Mg24$')
    plt.plot(x, h1, color='b', lw=1.5, label=r'$H1$')

    plt.xlim(0, x_lim)
    plt.xlabel(xlabel)
    plt.ylabel('mass fraction')
    plt.suptitle(title)
    plt.title('Abundance at model ' + str(model) + r' $M_i = $ ' +
              str(mass[0]) + r' $M_f = $ ' + str(mass[m_ind]) +
              r' $Teff = $ ' + str(10**teff[m_ind]))

    plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)

    if save is True:
        plt.savefig(name + '.pdf')
        return 'File save with name ' + name
    else:
        plt.show()


def get_mat_mcenter(path):
    '''
    Extract the initial and final mass, total star age, log L and Teff in the
    last model of evolution, and central masses of h1, he4, c12, n14, o16, ne20
    and mg24. Age in years

    Parameters
    ----------
    path : string
            LOGS directory path

    Returns
    -------
    total star age in years, log Teff, Teff
    '''

    fol = ms.MesaLogDir(log_path=path)
    h = fol.history
    age = h.data('star_age')
    logT = h.data('log_Teff')
    teff = 10**logT
    logL = h.data('log_L')
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

    return np.array([mi[0], mi[-1], age[-1], logL[-1],
                     teff[-1], c_h1[-1], c_he4[-1], c_c12[-1],
                     c_n14[-1], c_o16[-1], c_ne20[-1], c_mg24[0]])


def print_data(list_path):
    '''
    Print the initial and final mass, age, log L, Teff, center values of
    h1, he4, c12, n14, o16, ne20 and mg24

    Parameter
    --------
    list_path : list of strings
            list of all LOGS directories that are going to be analized
    '''

    np.set_printoptions(formatter={'float': '{: 0.5}'.format})

    print('Mi, Mf, age[yrs], log L, Teff, c_h1,' +
          'c_he4, c_c12, c_n14, c_o16, c_ne20, c_mg24')

    for i in list_path:
        print_data = get_mat_mcenter(i)
        print(print_data)


def plot_lum(folder, title='', save=False, name=None):
    '''
    Function to plot the Luminosity versus Age to check if occured
    cristalization

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
    age = hist.data('model_number')  # age = hist.data('star_age')
    mass = hist.data('star_mass')
    plt.figure(figsize=(10, 8))
    plt.plot(age, lum)
    plt.suptitle(r'$\log L$ X Model', fontsize=18)
    plt.title('Mi = ' + str(mass[0]) + ' - ' + title)
    plt.xlabel('model', fontsize=16)  # plt.xlabel('Age [years]', fontsize=16)
    plt.ylabel(r'$\log L$', fontsize=16)

    if save is True:
        plt.savefig(name + '.pdf')
        return 'File save with name ' + name
    else:
        plt.show()


def print_totalm(path):
    '''
    Print final mass and total mass of h1, he4, c12 and o16

    Parameter
    --------
    path : strings
            path for LOGS directory
    '''

    np.set_printoptions(formatter={'float': '{: 0.5}'.format})

    print('Total mass of isotopes (Msun)')
    print('mf,     h1,     he4,     c12,     o16,')
    fol = ms.MesaLogDir(path)
    hist = fol.history
    mass = hist.data('star_mass')
    t_h1 = hist.data('total_mass_h1')
    t_he = hist.data('total_mass_he4')
    t_c12 = hist.data('total_mass_c12')
    t_o16 = hist.data('total_mass_o16')
    print(mass[-1], t_h1[-1], t_he[-1], t_c12[-1], t_o16[-1])


def out_dir(out):
    '''
    check if the output path is relative or not and return the absolute path
    '''
    if out is None:
        return os.getcwd()
    elif os.path.isdir(out):
        return os.path.abspath(out)


def mod_abun(path, h_mass, name='abun.dat', output=None, profile_number=None):
    '''

    '''

    out = out_dir(output)
    path = ms.MesaLogDir(path)
    hist = path.history
    star_m = hist.data('star_mass')

    if profile_number is None:
        profiles = path.model_numbers
        model = profiles[-1]
        model_index = -1
    else:
        model = profile_number
        model_index = profile_number - 1
    prof = path.profile_data(model)

    headers = ['q', 'dm', 'xm', 'h1', 'he3', 'he4', 'c12',
               'c13', 'n13', 'n14', 'n15', 'o16', 'o17',
               'o18', 'f19', 'ne20', 'ne22', 'mg24', 'si28']

    data = np.array([prof.data('q'), prof.data('dm'), prof.data('xm'),
                     prof.data('h1'), prof.data('he3'), prof.data('he4'),
                     prof.data('c12'), prof.data('c13'), prof.data('n13'),
                     prof.data('n14'), prof.data('n15'), prof.data('o16'),
                     prof.data('o17'), prof.data('o18'), prof.data('f19'),
                     prof.data('ne20'), prof.data('ne22'), prof.data('mg24'),
                     prof.data('si28')])

    df = pd.DataFrame(data.T, columns=headers)
    m = h_mass * star_m[model_index]

    # Achar todos os indices em que he4 possui massa inferior a mx.
    # Para isto, primeiro eu faço a soma cumulativa da atmosfera até o nucleo
    # e identifico os valores menores que mx

    m_he4 = df['he4'].multiply(df['dm'], axis="index").multiply(
        5.027652086e-34).cumsum()

    # array com True onde m_he4 for  menor ou igual a m
    pos = m_he4.loc[:] <= m  # + 0.05 * m

    h1 = df['h1']
    he4 = df['he4']
    # Substituindo He por H
    h1[pos] = he4[pos]
    # Igualando He = 0 nos locais em que houve a substituicao
    he4[pos] = 0

    # Aplicando o filtro gaussiano
    gf_h1 = scipy.ndimage.gaussian_filter(h1, 3)
    gf_he4 = scipy.ndimage.gaussian_filter(he4, 3)

    new_data = df.copy()
    new_data['h1'] = gf_h1
    new_data['he4'] = gf_he4
    data_clean = new_data.drop(['q', 'dm'], axis=1)

    # Dados de colunas e linhas para o header do arquivo
    h1 = data_clean.shape[0]
    h2 = data_clean.shape[1] - 1
    header = str(h1) + ', ' + str(h2)

    # Salvando os dados
    np.savetxt(os.path.join(out, str(name)), data_clean, fmt='%1.15f',
               header=header, comments='')

    # Print da quantidade de H1 adicionada
    print(np.sum(data_clean.h1 * df['dm'] * 5.027652086e-34))
