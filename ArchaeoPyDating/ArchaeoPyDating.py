import numpy as np
from numpy import pi, cos, sin, tan, arctan, arccos, sqrt, exp, radians, degrees
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, griddata
import shapefile
import warnings
from .reference_curves import *
from past.utils import old_div
import os
from pyshtools import legendre

warnings.filterwarnings("ignore")

# directory of this script
dir_path = os.path.dirname(os.path.realpath(__file__))


class Model:
    def __init__(self, name):
        self.name = name
        self.coeff = np.loadtxt(os.path.join(dir_path, 'curves', 'gmodel', global_models[name][0]))
        self.ecoeff = np.loadtxt(os.path.join(dir_path, 'curves', 'gmodel', global_models[name][1]))
        self.t = np.loadtxt(os.path.join(dir_path, 'curves', 'gmodel', global_models[name][2]))

    def normal_distribute(self):
        '''Return random Gauss coefficients acording to a normal distribution'''
        return np.random.normal(self.coeff, self.ecoeff)

    def psvc(self, lat, lon):
        '''Create a psvc at lat, lon.'''

        nmax = int((-2 + np.sqrt(4 + 4 * self.coeff.shape[1])) / 2)
        coord = np.array([[np.radians(90 - lat), np.radians(lon), 6371.2]])  # Make a 2D array
        A = frechet_basis(coord, nmax)
        A = A.squeeze(0) # Remove the first dimension
        
        XYZ = np.dot(A, self.coeff.T)
        X, Y, Z = XYZ

        D, I = cart2dir(X, Y, Z)
        F = np.sqrt(X ** 2 + Y ** 2 + Z ** 2)
        H = np.sqrt(X ** 2 + Y ** 2)

        eXYZ = np.sqrt( np.dot( A**2,  self.ecoeff.T**2  ) )
        eX, eY, eZ = eXYZ

        eD = np.degrees(  (1/H**2)*np.sqrt( (Y*eX)**2 + (X*eY)**2 )   )
        eI = np.degrees(  (1/(F**2*H))*np.sqrt( (X*Z*eX)**2 + (Y*Z*eY)**2 + (H**2*eZ)**2 ) )
        eF = (1/F)*np.sqrt( (X*eX)**2 + (Y*eY)**2 + (Z*eZ)**2 )

        return D, I, F, eD, eI, eF

    def field_t(self, t, lat, lon):
        D, I, F, eD, eI, eF = self.psvc(lat, lon)
        matriz = np.column_stack([D, I, F])
        cmatriz = CubicSpline(self.t, matriz)
        Dt, It, Ft = cmatriz(t)
        return Dt, It, Ft


class RegionalModel:
    def __init__(self, name):
        self.name = name
        self.matriz = np.loadtxt(os.path.join(dir_path, 'curves', 'rmodel', regional_models[name]))
        self.or_lat, self.or_lon = self.matriz[0, 0], self.matriz[0, 1]
        self.angular_distance = self.matriz[0, 2]
        self.matriz = np.delete(self.matriz, 0, 0)  # para quitar la cabecera
        self.t = np.unique(self.matriz[:, 0])

    def psvc(self, nlat, nlon):
        Dc = np.zeros(len(self.t));
        Ic = Dc.copy();
        Fc = Dc.copy();
        eDc = Dc.copy();
        eIc = Dc.copy();
        eFc = Dc.copy();
        Xc = [Dc, Ic, Fc];
        eXc = [eDc, eIc, eFc]
        for i, ti in enumerate(self.t):
            y = self.matriz[self.matriz[:, 0] == ti, :]
            y = y[distance(y[:, 1], y[:, 2], nlat, nlon) <= 3,
                :]  # para hacer el calculo rapido elimino aquello + lejos de 3º
            lat, lon = y[:, 1], y[:, 2]
            D, I, F = y[:, 3], y[:, 4], y[:, 7]
            eD, eI, eF = y[:, 5], y[:, 6], y[:, 8]
            X = [D, I, F]
            eX = [eD, eI, eF]

            for x, ex, xc, exc in zip(X, eX, Xc, eXc):
                nx = griddata((lat, lon), x, (nlat, nlon), method='cubic')
                nex = griddata((lat, lon), ex, (nlat, nlon), method='cubic')
                xc[i] = nx
                exc[i] = nex
        return Dc, Ic, Fc, eDc / 1.96, eIc / 1.96, eFc / 1.96


class Data:
    def __init__(self, D=0, I=0, a95=0, F=0, eF=0, lat=None, lon=None, sitename='Site'):
        self.D = D
        self.I = I
        self.a95 = a95
        self.eI = self.a95 / 2.45  # (65%) Lanos et al. 2005; Suttie and Nilsson, 2019; Nilsson and Suttie, 2021
        self.eD = self.eI / cos(radians(self.I))
        self.F = F
        self.eF = eF
        self.lat = lat
        self.lon = lon
        self.sitename = sitename

    def cvp(self, nlat, nlon):
        # radians
        lat, lon, D, I, nlat, nlon = radians([self.lat, self.lon, self.D, self.I, nlat, nlon])
        self.Dold, self.Iold = np.degrees(D), np.degrees(I)
        # from lat 2 colat
        clat, nclat = pi / 2 - np.array([lat, nlat])
        # dipolar coef
        g10 = - (cos(D) * cos(I) * sin(clat) + 0.5 * sin(I) * cos(clat))
        g11 = cos(D) * cos(I) * cos(clat) * cos(lon) + sin(D) * cos(I) * sin(lon) - 0.5 * sin(I) * sin(clat) * cos(lon)
        h11 = - (-cos(D) * cos(I) * cos(clat) * sin(lon) + sin(D) * cos(I) * cos(lon) + 0.5 * sin(I) * sin(clat) * sin(
            lon))
        # comp. x,y,z at the new point
        x = -g10 * sin(nclat) + cos(nclat) * (g11 * cos(nlon) + h11 * sin(nlon))
        y = g11 * sin(nlon) - h11 * cos(nlon)
        z = -2 * (g10 * cos(nclat) + sin(nclat) * (g11 * cos(nlon) + h11 * sin(nlon)))
        # news D, I
        Dn, In = cart2dir(x, y, z)
        self.D = Dn
        self.I = In

    def vadm(self, nlat):
        # Translate the intensity data using the Virtual Axial Dipole Moment
        # radians
        lat, nlat = radians([self.lat, nlat])
        # from lat 2 colat
        clat, nclat = pi / 2 - np.array([lat, nlat])
        Fn = self.F * sqrt((1 + 3 * cos(nclat) ** 2) / (1 + 3 * cos(clat) ** 2))
        self.Fold = self.F
        self.F = Fn

    def vdm(self):
        # Translate the intensity data using the Virtual Dipole Moment
        # you know inclination and intensity in one point
        # and only inclination in a new point
        # radians
        I, In = radians(np.array([self.Iold, self.I]))
        clat, nclat = arctan(2 / tan([I, In]))
        Fn = self.F * sqrt((1 + 3 * cos(nclat) ** 2) / (1 + 3 * cos(clat) ** 2))
        self.Fold = self.F
        self.F = Fn


class Dating:
    def __init__(self, dato, curva):
        self.dato = dato
        self.curva = curva

    def datingX(self, indicador):  # indicator: 'D', 'I', 'F'
        if indicador == 'D':
            inter = np.arange(-pi, pi, radians(0.1))
            X, eX = self.dato.D, self.dato.eD
            Xc, eXc = self.curva.D, self.curva.eD
            X, eX, Xc, eXc = radians(X), radians(eX), radians(Xc), radians(eXc)
        if indicador == 'I':
            inter = np.arange(-pi / 2, pi / 2, radians(0.1))
            X, eX = self.dato.I, self.dato.eI
            Xc, eXc = self.curva.I, self.curva.eI
            X, eX, Xc, eXc = radians(X), radians(eX), radians(Xc), radians(eXc)
        if indicador == 'F':
            inter = np.arange(0, 100, 0.1)
            X, eX = self.dato.F, self.dato.eF
            Xc, eXc = self.curva.F, self.curva.eF

        # data fun. probability
        fd1 = exp(- (inter - X) ** 2 / (2 * eX ** 2)) / (eX * sqrt(2 * pi))
        # curve fun. probability (depends on time)
        # product data fun. and curve fun.
        area = []
        for xc, exc in zip(Xc, eXc):
            fd2 = exp(- (inter - xc) ** 2 / (2 * exc ** 2)) / (exc * sqrt(2 * pi))
            pro = fd1 * fd2
            areai = np.trapz(pro, inter)  # trapezoidal integration
            area = np.append(area, areai)
        # normalization
        area = area / np.trapz(area, self.curva.t)
        # interpolation 1yr
        t = np.arange(self.curva.t[0], self.curva.t[-1] + 1, 1)
        cs = CubicSpline(self.curva.t, abs(area))
        z = cs(t)
        return z

    def zcomb(self, *args):
        aux = 1
        for z in args:
            aux = aux * z
        t = np.arange(self.curva.t[0], self.curva.t[-1] + 1, 1)
        z = aux / np.trapz(aux, t)
        return z

    def pb(self, z, gp):
        '''z from Dating.datingX()'''
        H = np.linspace(min(z), max(z), 501)
        t = np.arange(self.curva.t[0], self.curva.t[-1] + 1, 1)
        area = []
        for h in H:
            areai = np.trapz(z * [z >= h], t)
            area = np.append(area, areai)
        pm = H[area <= gp / 100]
        pM = H[area >= gp / 100]
        if (pm.size == 0) | (pM.size == 0):
            print("An archaeomagnetic dating with these data is not possible")
            return
        h = (pm[0] + pM[-1]) / 2
        return h

    def pb_h(self, z, h):
        '''z from Dating.datingX() and h from Dating.pb()'''
        # intersection of z and h
        # new function and its solutions
        t = np.arange(self.curva.t[0], self.curva.t[-1] + 1, 1)
        f = CubicSpline(t, z - h)
        sol = f.roots()
        # no solutions outside the interval
        sol = sol[(sol > t[0]) & (sol < t[-1])]
        # search if difference of sol are positive or negative
        # to know if they represent lower or upper bounds
        if f(sol[0], 1) > 0:
            tmin = sol[::2]
            tmax = sol[1::2]
            # if the last interval does not end, I end it at the time maximum.
            if len(tmin) > len(tmax):
                tmax = np.append(tmax, t[-1])
        if f(sol[0], 1) < 0:
            tmin = sol[1::2]
            tmax = sol[::2]
            # if the first interval does not start, I start it at the minimum time interval
            if len(tmin) < len(tmax):
                tmin = np.append(t[0], tmin)
        # it could be that len of tmin == tmax because it neither begins nor ends.
        if tmin[0] > tmax[0]:
            tmax = np.append(tmax, t[-1])
            tmin = np.append(t[0], tmin)
        # datings
        dat = {}
        for i in range(len(tmin)):
            dat["Dating result {0}".format(i + 1)] = [int(round(tmin[i] / 5) * 5), int(round(tmax[i] / 5) * 5)]
        result = {key: [dat[key][0], dat[key][1]] for key in dat}
        return result

    def plot(self, zD=None, zI=None, zF=None, z=None,
        hD=None, hI=None, hF=None, h=None):
    
        fig = plt.figure(figsize=(18, 12), dpi=600)
        plt.rcParams.update({'font.size': 14})
        # create grid with all
        gs = fig.add_gridspec(nrows=3, ncols=2)
        # within that grid create three grids for the left side.
        gsD = gs[0, 0].subgridspec(2, 1, hspace=.0, height_ratios=[3, 1.5])
        gsI = gs[1, 0].subgridspec(2, 1, hspace=0, height_ratios=[3, 1.5])
        gsF = gs[2, 0].subgridspec(2, 1, hspace=0, height_ratios=[3, 1.5])
        # within each of these grids, put two ax
        axX = [None, None, None]
        axfX = [None, None, None]
        gsX = [gsD, gsI, gsF]
        for i, (gsXi) in enumerate(gsX):
            axX[i] = fig.add_subplot(gsXi[0])
            plt.tick_params('x', labelbottom=False, direction='in')
            axfX[i] = fig.add_subplot(gsXi[1], sharex=axX[i])

        # ax of f comb
        axfcomb = fig.add_subplot(gs[0, 1], sharex=axfX[0])

        # grid for the map
        gs_info = gs[1, 1].subgridspec(1, 1)
        axmapa = fig.add_subplot(gs_info[0])

        # ax for dating
        axdat = fig.add_subplot(gs[2, 1])

        # define the things to draw
        ylabel = ['Declination [$\degree$]', 'Inclination [$\degree$]', 'Intensity [$\mu$T]']
        X = [self.dato.D, self.dato.I, self.dato.F]
        eX = [self.dato.eD, self.dato.eI, self.dato.eF]
        Xc = [self.curva.D, self.curva.I, self.curva.F]
        eXc = [self.curva.eD, self.curva.eI, self.curva.eF]
        zX = [zD, zI, zF]
        tz = np.arange(self.curva.t[0], self.curva.t[-1] + 1, 1)
        H = [hD, hI, hF]

        # plot the curves with the data
        for ax, zx, x, ex, xc, exc, label in zip(axX, zX, X, eX, Xc, eXc, ylabel):
            if zx is not None:
                ax.plot(self.curva.t, xc, 'r-')
                ax.fill_between(self.curva.t, (xc + 2 * exc), (xc - 2 * exc), color='r', alpha=0.1)
                ax.axhline(x, color='g')
                ax.axhspan(x - 2 * ex, x + 2 * ex, color='g', alpha=0.1)
                ax.set_xlim([self.curva.t[0], self.curva.t[-1]])
                ax.set_ylabel(label)
            if zx is None:
                ax.remove()

        # plot the PDFs
        for ax, zx, hx in zip(axfX, zX, H):
            if zx is not None:
                ax.plot(tz, zx, 'k-')
                ax.set_ylim(bottom=0)
                ax.set_ylabel('PDF', multialignment='center')
                ax.set_xlabel('Year')
                ax.axhline(hx)
                ax.fill_between(tz, zx, where=zx > hx, color='gold')
                ax.set_yticklabels([])
            if zx is None:
                ax.remove()

        # plot the combined PDF
        axfcomb.plot(tz, z, 'k-')
        axfcomb.axhline(h)
        axfcomb.fill_between(tz, z, where=z > h, color='gold')
        axfcomb.set_ylim(bottom=0)
        axfcomb.set_ylabel('Combined PDF', multialignment='center')
        axfcomb.set_xlabel('Year')

        sf = shapefile.Reader('external/ne_50m_coastline.shp')
        for shape in sf.shapes():
            points = shape.points 
            x, y = zip(*points)  
            axmapa.plot(x, y, 'k-', linewidth=1) 
        # Add points for the locations
        axmapa.scatter(self.dato.lon, self.dato.lat, color='blue', s=100, label='Data', zorder=5)
        axmapa.scatter(self.curva.lon, self.curva.lat, color='orange', s=100, marker='*', label='Curve', zorder=5)
#
        axmapa.set_xlim(self.dato.lon - 30, self.dato.lon + 30)
        axmapa.set_ylim(self.dato.lat - 20, self.dato.lat + 20)
        axmapa.legend()

        # hide ax with text
        axdat.axis('off')
        
        
        plt.tight_layout()
        plt.show()
        return fig


class Curve:
    def __init__(self, local=None, gmodel=None, rmodel=None, lat=None, lon=None, newpsvc=None):
        if local:
            self.name = local
            self.matriz = np.loadtxt(os.path.join(dir_path, 'curves', 'local', local_curves[local][2]))
            self.lat, self.lon = local_curves[local][0:2]
            if self.matriz.shape[1] == 7:
                self.t, self.D, self.I, self.eD, self.eI, self.F, self.eF = self.matriz.T
            if self.matriz.shape[1] == 5:
                self.t, self.D, self.I, self.eD, self.eI = self.matriz.T
                self.F, self.eF = np.zeros(self.matriz.shape[0]), np.zeros(self.matriz.shape[0])
            if self.matriz.shape[1] == 3:
                self.t, self.F, self.eF = self.matriz.T
                self.D, self.I, self.eD, self.eI = np.zeros(self.matriz.shape[0]), np.zeros(self.matriz.shape[0]), np.zeros(self.matriz.shape[0]), np.zeros(self.matriz.shape[0])

        if gmodel:
            mod = Model(gmodel)
            self.name = gmodel + '  (global model)'
            self.t = mod.t
            self.D, self.I, self.F, self.eD, self.eI, self.eF = mod.psvc(lat, lon)
            self.lat, self.lon = lat, lon
            self.matriz = np.column_stack([self.t, self.D, self.I, self.eD, self.eI, self.F, self.eF])

        if rmodel:
            mod = RegionalModel(rmodel)
            self.name = rmodel + '  (regional model)'
            self.t = mod.t
            self.D, self.I, self.F, self.eD, self.eI, self.eF = mod.psvc(lat, lon)
            self.lat, self.lon = lat, lon
            self.matriz = np.column_stack([self.t, self.D, self.I, self.eD, self.eI, self.F, self.eF])
            self.or_lat, self.or_lon = mod.or_lat, mod.or_lon
            self.angular_distance = mod.angular_distance

        if newpsvc:

            aux = np.asarray(newpsvc)
            aux2 = np.zeros([len(aux), len(np.fromstring(aux[0], dtype=float, sep=' '))])
            for i in range(len(aux)):
                aux2[i] = np.fromstring(aux[i], dtype=float, sep=' ')
            self.matriz = aux2

            self.name = 'Private curve'
            self.lat, self.lon = lat, lon
            if self.matriz.shape[1] == 7:
                self.t, self.D, self.I, self.eD, self.eI, self.F, self.eF = self.matriz.T
            if self.matriz.shape[1] == 5:
                self.t, self.D, self.I, self.eD, self.eI = self.matriz.T
            if self.matriz.shape[1] == 3:
                self.t, self.F, self.eF = self.matriz.T

        # interpolte to 10 yr
        t = np.arange(self.t[0], self.t[-1] + 10, 10)
        d = CubicSpline(self.t, self.D);
        i = CubicSpline(self.t, self.I);
        f = CubicSpline(self.t, self.F)
        ed = CubicSpline(self.t, self.eD);
        ei = CubicSpline(self.t, self.eI);
        ef = CubicSpline(self.t, self.eF)
        self.D, self.I, self.F = d(t), i(t), f(t)
        self.eD, self.eI, self.eF = ed(t), ei(t), ef(t)
        self.t = t

    def int_temp(self, tmin, tmax):  # if a time interval is defined
        self.D = self.D[(self.t >= tmin) & (self.t <= tmax)];
        self.eD = self.eD[(self.t >= tmin) & (self.t <= tmax)]
        self.I = self.I[(self.t >= tmin) & (self.t <= tmax)];
        self.eI = self.eI[(self.t >= tmin) & (self.t <= tmax)]
        self.F = self.F[(self.t >= tmin) & (self.t <= tmax)];
        self.eF = self.eF[(self.t >= tmin) & (self.t <= tmax)]
        self.t = self.t[(self.t >= tmin) & (self.t <= tmax)]


def distance(lat1, lon1, lat2, lon2):
    '''Distance between two points on sphere in degrees'''
    lat1, lon1, lat2, lon2 = radians(lat1), radians(lon1), radians(lat2), radians(lon2)
    angular_distance = arccos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2))
    return degrees(angular_distance)


def cart2dir(X, Y, Z):
    D = np.degrees(np.arctan(Y / X))
    I = np.degrees(np.arctan(Z / np.sqrt(X ** 2 + Y ** 2)))
    return D, I

def frechet_basis(loc: np.ndarray,
                  maxdegree: int
                  ) -> np.ndarray:
    """
    Calculates the frechet matrix for the given stations and maximum degree
    Parameters
    ----------
    loc
        coordinates of stations. Each row contains:
            colatitude in radians
            longitude in radians
            radius in km
    maxdegree
        Maximum spherical degree

    Returns
    -------
    frechxyz
        size= stations X 3 X nr_coeffs matrix: contains the frechet coefficients
            first dx, then dy, then dz
    """
    schmidt_total = int((maxdegree+1) * (maxdegree+2) / 2)
    frechxyz = np.zeros(((maxdegree+1)**2 - 1, 3, len(loc)))
    schmidt_p = np.zeros((len(loc), schmidt_total))
    schmidt_dp = np.zeros((len(loc), schmidt_total))
    for i, coord in enumerate(loc):
        schmidt_p[i], schmidt_dp[i] = \
            legendre.PlmSchmidt_d1(maxdegree, np.cos(coord[0]))
        schmidt_dp[i] *= -np.sin(coord[0])
    counter = 0
    # dx, dy, dz in separate rows to increase speed
    for n in range(1, maxdegree+1):
        index = int(n * (n+1) / 2)
        mult_factor = (6371.2 / loc[:, 2]) ** (n+1)
        # first g_n^0
        frechxyz[counter, 0] = mult_factor * schmidt_dp[:, index]
        # frechxyz[counter, 1] = 0
        frechxyz[counter, 2] = -mult_factor * (n+1) * schmidt_p[:, index]
        counter += 1
        for m in range(1, n+1):
            # Then the g-elements
            frechxyz[counter, 0] = mult_factor\
                * np.cos(m * loc[:, 1]) * schmidt_dp[:, index+m]
            frechxyz[counter, 1] = m / np.sin(loc[:, 0]) * mult_factor\
                * np.sin(m * loc[:, 1]) * schmidt_p[:, index+m]
            frechxyz[counter, 2] = -mult_factor * (n+1)\
                * np.cos(m * loc[:, 1]) * schmidt_p[:, index+m]
            counter += 1
            # Now the h-elements
            frechxyz[counter, 0] = mult_factor\
                * np.sin(m * loc[:, 1]) * schmidt_dp[:, index+m]
            frechxyz[counter, 1] = -m / np.sin(loc[:, 0]) * mult_factor\
                * np.cos(m * loc[:, 1]) * schmidt_p[:, index+m]
            frechxyz[counter, 2] = -mult_factor * (n+1)\
                * np.sin(m * loc[:, 1]) * schmidt_p[:, index+m]
            counter += 1
    # transpose to get frechet matrix
    return np.swapaxes(frechxyz, 0, 2)
