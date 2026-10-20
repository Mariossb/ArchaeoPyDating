import numpy as np
from numpy import pi, cos, sin, tan, arctan, arccos, sqrt, exp, radians, degrees
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, griddata
from scipy.special import ive
import shapefile
from .reference_curves import *
import os
from functools import lru_cache

__all__ = ['Curve', 'Data', 'Dating', 'Model', 'RegionalModel', 'available_psvc']

# directory of this script
dir_path = os.path.dirname(os.path.realpath(__file__))


class Model:
    def __init__(self, name):
        self.name = name
        self.coeff = np.loadtxt(os.path.join(dir_path, 'curves', 'gmodel', global_models[name][0]))
        self.ecoeff = np.loadtxt(os.path.join(dir_path, 'curves', 'gmodel', global_models[name][1]))
        self.t = np.loadtxt(os.path.join(dir_path, 'curves', 'gmodel', global_models[name][2]))

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

        eD, eI, eF = np.zeros(len(X)), np.zeros(len(X)), np.zeros(len(X))
        for k in range(len(self.t)):
            # covariance of X, Y, Z: A Sg A^T, with Sg the covariance of the
            # coefficients. With only variances, Sg is diagonal and
            # A Sg A^T = (A * variances) A^T
            S = (A * self.ecoeff[k] ** 2) @ A.T

            # jacobian of (D, I, F) with respect to (X, Y, Z); D and I in radians
            x, y, z, f, h = X[k], Y[k], Z[k], F[k], H[k]
            J = np.array([[-y / h**2,            x / h**2,            0       ],
                          [-x * z / (h * f**2),  -y * z / (h * f**2), h / f**2],
                          [x / f,                y / f,               z / f   ]])

            # covariance of (D, I, F)
            C = J @ S @ J.T
            eD[k] = np.degrees(np.sqrt(C[0, 0]))
            eI[k] = np.degrees(np.sqrt(C[1, 1]))
            eF[k] = np.sqrt(C[2, 2])
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

    def psvc(self, lat, lon):
        # only take grid points within 3 degrees of the site
        near = distance(self.matriz[:,1], self.matriz[:,2], lat, lon) <= 3
        sub = self.matriz[near]
        if len(sub) == 0:
            raise ValueError(f"The site ({lat}, {lon}) is outside the area "
            f"of the regional model {self.name}")
        
        # columns of the file: D, I, F, eD, eI, eF
        columns = [3, 4, 7, 5, 6, 8]
        result = np.zeros((len(self.t), 6))
        for i, ti in enumerate(self.t):
            y = sub[sub[:, 0] == ti]
            result[i] = griddata((y[:, 1], y[:, 2]), y[:, columns], (lat, lon), method='cubic')
        if np.isnan(result).any():       
            raise ValueError(f"The site ({lat}, {lon}) is outside the area "
            f"of the regional model {self.name}")

        D, I, F, eD, eI, eF = result.T
        return D, I, F, eD / 1.96, eI / 1.96, eF / 1.96


class Data:
    def __init__(self, D=None, I=None, a95=None, F=None, eF=None, lat=None, lon=None, sitename='Site'):
        self.D = D
        self.I = I
        self.a95 = a95
        self.F = F
        self.eF = eF
        self.lat = lat
        self.lon = lon
        self.sitename = sitename

        # direction errors (65%) Lanos et al. 2005; Suttie and Nilsson, 2019; Nilsson and Suttie, 2021
        self.eI = a95 / 2.45  if a95 is not None else None
        if self.eI is not None and I is not None:
            self.eD = self.eI / cos(radians(I))
        else:
            self.eD = None   
        

    def cvp(self, nlat, nlon):
        if self.D is None or self.I is None:
            raise ValueError("cvp needs both declination and inclination")
        self.Dold, self.Iold = self.D, self.I
        self.D, self.I = cvp_direction(self.D, self.I, self.lat, self.lon, nlat, nlon)
        if self.eI is not None:
            self.eD = self.eI / cos(radians(self.I))  # same a95 cone, new inclination

    def vadm(self, nlat):
        # Translate the intensity data using the Virtual Axial Dipole Moment
        if self.F is None:
            raise ValueError("vadm needs intensity")
        self.Fold = self.F
        self.F = vdm_intensity(self.F, self.Iold, self.I)

    def vdm(self):
        # Translate the intensity data using the Virtual Dipole Moment
        # you know inclination and intensity in one point
        # and only inclination in a new point
        if self.F is None or self.I is None or not hasattr(self, 'Iold'):
            raise ValueError("vdm needs intensity and inclination, relocated first with cvp")
        self.Fold = self.F
        self.F = vadm_intensity(self.F, self.lat, nlat)


class Dating:
    def __init__(self, dato, curva):
        self.dato = dato
        self.curva = curva
        self.t = self.curva.t

    def datingX(self, indicador):  # indicator: 'D', 'I', 'F'
        if indicador == 'D':
            X, eX = self.dato.D, self.dato.eD
            Xc, eXc = self.curva.D, self.curva.eD
        if indicador == 'I':
            X, eX = self.dato.I, self.dato.eI
            Xc, eXc = self.curva.I, self.curva.eI
        if indicador == 'F':
            X, eX = self.dato.F, self.dato.eF
            Xc, eXc = self.curva.F, self.curva.eF
        
        if X is None or eX is None:
            raise ValueError(f"The data has no value or error for {indicador}")
        if Xc is None or eXc is None:
            raise ValueError(f"The curve has no {indicador}")
        if indicador == 'D':
            Xc = X - ((X - Xc + 180) % 360 - 180)  # curve values on the same turn as X [0,360] or [-180,180]
        
        diff = X - Xc
        var = eX**2 + eXc**2
        z = exp(-diff ** 2 / (2 * var)) / sqrt(2 * pi * var)
        z = z / np.trapezoid(z, self.t)
        return z

    def datingDI(self):
        d, c = self.dato, self.curva
        if d.D is None or d.I is None or d.a95 is None:
            raise ValueError("The data has no declination, inclination or a95")
        if c.D is None or c.I is None:
            raise ValueError("The curve has no declination and inclination")

        kd = (140/d.a95)**2  #a95 in degrees, kd in (rad^-1)^2
        kc = 2/( np.radians(c.eI)**2 + (np.radians(c.eD) * np.cos(np.radians(c.I))) ** 2) 

        Id, Ic = np.radians(d.I), np.radians(c.I)
        cos_theta = np.sin(Id) * np.sin(Ic) + np.cos(Id) * np.cos(Ic) * np.cos(np.radians(d.D - c.D))
        R = np.sqrt(kd ** 2 + kc ** 2 + 2 * kd * kc * cos_theta)

        logL = np.log(kc) - log_sinh(kc) + log_sinh(R) - np.log(R)
        z = np.exp(logL - logL.max())
        return z / np.trapezoid(z, self.t)

    def datingI(self):
        '''PDF of an inclination-only data with the marginal of the Fisher
        distribution over the unknown declination (McFadden and Reid, 1982)'''
        d, c = self.dato, self.curva
        if d.I is None or d.a95 is None:
            raise ValueError("The data has no inclination or a95")
        if c.I is None:
            raise ValueError("The curve has no inclination")

        kd = (140 / d.a95) ** 2 #a95 in degrees, kd in (rad^-1)^2
        kc = 2/( np.radians(c.eI)**2 + (np.radians(c.eD) * np.cos(np.radians(c.I))) ** 2) 
        k = 1 / (1 / kd + 1 / kc)

        theta = radians(90 - d.I)
        theta0 = radians(90 - c.I)
        x = k * sin(theta0) * sin(theta)

        # log of f; sin(theta) and the factor 2 are the same in every epoch.
        # ive(0, x) = I0(x) exp(-x), so log I0(x) = log(ive(0, x)) + x
        logL = np.log(k) - log_sinh(k) + k * cos(theta0) * cos(theta) + np.log(ive(0, x)) + x
        z = exp(logL - logL.max())
        return z / np.trapezoid(z, self.t)

    def zcomb(self, *args):
        z = np.ones_like(self.t)
        for zi in args:
            z *= zi
        # normalizar
        z = z / np.trapezoid(z, self.t)
        return z

    def pb(self, z, gp):
        '''Height h of the PDF whose region z>=h encloses gp% of the probability'''
        if not np.all(np.isfinite(z)) or np.max(z)<=0:
            raise ValueError("Dating not possible with these data")
        
        order = np.argsort(z)[::-1] #indices of z form highest to lowest
        cum = np.cumsum(z[order]) #probability accumulated in that order
        cum /= cum[-1] #from 0 to 1
        k = np.searchsorted(cum, gp/100)
        h = z[order[k]]
        return h

    def pb_h(self, z, h, step=5):
        mask = z >= h

        intervals = []
        start = None

        for i in range(len(self.t)):
            if mask[i] and start is None:
                start = self.t[i]
            elif not mask[i] and start is not None:
                end = self.t[i-1]
                intervals.append((start, end))
                start = None

        # close last interval if open
        if start is not None:
            intervals.append((start, self.t[-1]))

        # round results and join close intervals
        rounded = []
        for t1, t2 in intervals:
            t1 = step * np.floor(t1/step)
            t2 = step * np.ceil(t2/step)
            if rounded and t1 <= rounded[-1][1]:
                rounded[-1][1] = t2
            else:
                rounded.append([t1, t2])

        result = {}
        for i, (t1, t2) in enumerate(rounded):
            result[f"Dating result {i+1}"] = [int(t1), int(t2)]

        return result

    @plt.rc_context({'font.size': 14})
    def plot(self, zD=None, zI=None, zF=None, z=None,
        hD=None, hI=None, hF=None, h=None):
    
        fig = plt.figure(figsize=(18, 12))
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
            axX[i].tick_params('x', labelbottom=False, direction='in')
            axfX[i] = fig.add_subplot(gsXi[1], sharex=axX[i])

        # ax of f comb
        axfcomb = fig.add_subplot(gs[0, 1], sharex=axfX[0])

        # grid for the map
        gs_info = gs[1, 1].subgridspec(1, 1)
        axmapa = fig.add_subplot(gs_info[0])

        # ax for dating
        axdat = fig.add_subplot(gs[2, 1])

        # define the things to draw
        ylabel = [r'Declination [$\degree$]', r'Inclination [$\degree$]', r'Intensity [$\mu$T]']
        X = [self.dato.D, self.dato.I, self.dato.F]
        eX = [self.dato.eD, self.dato.eI, self.dato.eF]
        Xc = [self.curva.D, self.curva.I, self.curva.F]
        eXc = [self.curva.eD, self.curva.eI, self.curva.eF]
        zX = [zD, zI, zF]
        tz = self.t
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

        coastlon, coastlat = coastlines()
        axmapa.plot(coastlon, coastlat, 'k-', linewidth=1)

        # Add points for the locations
        axmapa.scatter(self.dato.lon, self.dato.lat, color='blue', s=100, label='Data', zorder=5)
        axmapa.scatter(self.curva.lon, self.curva.lat, color='orange', s=100, marker='*', label='Curve', zorder=5)

        axmapa.set_xlim(self.dato.lon - 30, self.dato.lon + 30)
        axmapa.set_ylim(self.dato.lat - 20, self.dato.lat + 20)
        axmapa.legend()

        # hide ax with text
        axdat.axis('off')

        result = self.pb_h(z, h)

        if len(result)<=6:

            text = '\n\n'+'Dating results:'+'\n\n'
            counter = 0
            for key, value in result.items():
                m = BCADyear(value[0])
                M = BCADyear(value[1])
                text += m+' - '+M+';  '
                counter += 1
                if counter % 3 ==0: text += '\n'

            axdat.text(0.5, 0.2, text, fontsize=20, fontweight='bold', ha='center')
        
        
        fig.tight_layout()
        return fig


class Curve:
    def _read_columns(self, source):
        '''Assign the columns of a curve: 7 (t, D, I, eD, eI, F, eF),
        5 (t, D, I, eD, eI) or 3 (t, F, eF). Missing elements are None.'''
        matriz = np.atleast_2d(np.loadtxt(source))
        ncol = matriz.shape[1]
        self.D = self.I = self.eD = self.eI = self.F = self.eF = None
        if ncol == 7:
            self.t, self.D, self.I, self.eD, self.eI, self.F, self.eF = matriz.T
        elif ncol == 5:
            self.t, self.D, self.I, self.eD, self.eI = matriz.T
        elif ncol == 3:
            self.t, self.F, self.eF = matriz.T
        else:
            raise ValueError(f"A curve must have 7, 5 or 3 columns, not {ncol}")
        self.matriz = matriz

    def __init__(self, local=None, gmodel=None, rmodel=None, lat=None, lon=None, newpsvc=None):
        if sum(x is not None for x in (local, gmodel, rmodel, newpsvc)) != 1:
            raise ValueError("Give one curve/model")
        if local:
            if not local_curves.get(local):
                raise ValueError(f"Unknown local curve: {local}")
            self.name = local
            self.lat, self.lon = local_curves[local][0:2]
            self._read_columns(os.path.join(dir_path, 'curves', 'local', local_curves[local][2]))
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
            lines = [l.decode() if isinstance(l, bytes) else l for l in newpsvc]
            self.name = 'Private curve'
            self.lat, self.lon = lat, lon
            self._read_columns(lines)

        # interpolate to 1 yr
        t = np.arange(self.t[0], self.t[-1] + 0.5, 1)
        for name in ('D', 'I', 'F', 'eD', 'eI', 'eF'):
            values = getattr(self, name)
            if values is not None:
                setattr(self, name, CubicSpline(self.t, values)(t))
        self.t = t

    def int_temp(self, tmin, tmax):  # if a time interval is defined
        mask = (self.t >= tmin) & (self.t <= tmax)
        for name in ('D', 'I', 'F', 'eD', 'eI', 'eF'):
            values = getattr(self, name)
            if values is not None:
                setattr(self, name, values[mask])
        self.t = self.t[mask]


def distance(lat1, lon1, lat2, lon2):
    '''Distance between two points on sphere in degrees'''
    lat1, lon1, lat2, lon2 = radians(lat1), radians(lon1), radians(lat2), radians(lon2)
    angular_distance = arccos(np.clip(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2), -1, 1))
    return degrees(angular_distance)


def cart2dir(X, Y, Z):
    D = np.degrees(np.arctan2(Y, X))
    I = np.degrees(np.arctan2(Z, np.sqrt(X ** 2 + Y ** 2)))
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
        schmidt_p[i], schmidt_dp[i] = schmidt_legendre(maxdegree, coord[0])
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

def schmidt_legendre(nmax, theta):
    '''Schmidt semi-normalized associated Legendre functions P_n^m(cos theta)
    and their derivatives with respect to theta, for 0 <= m <= n <= nmax.
    Returned packed as in pyshtools: index n(n+1)/2 + m.'''
    c, s = np.cos(theta), np.sin(theta)
    P = np.zeros((nmax + 1, nmax + 1))
    dP = np.zeros((nmax + 1, nmax + 1))
    P[0, 0] = 1.0
    for n in range(1, nmax + 1):
        # diagonal term P_n^n from P_{n-1}^{n-1}
        k = 1.0 if n == 1 else np.sqrt((2 * n - 1) / (2 * n))
        P[n, n] = k * s * P[n - 1, n - 1]
        dP[n, n] = k * (c * P[n - 1, n - 1] + s * dP[n - 1, n - 1])
        # other orders m < n from degrees n-1 and n-2
        for m in range(n):
            a = (2 * n - 1) / np.sqrt(n * n - m * m)
            b = np.sqrt(((n - 1) ** 2 - m * m) / (n * n - m * m))
            P[n, m] = a * c * P[n - 1, m] - (b * P[n - 2, m] if n > 1 else 0.0)
            dP[n, m] = a * (c * dP[n - 1, m] - s * P[n - 1, m]) - (b * dP[n - 2, m] if n > 1 else 0.0)
    rows, cols = np.tril_indices(nmax + 1)
    return P[rows, cols], dP[rows, cols]

@lru_cache(maxsize=1)
def coastlines():
    sf = shapefile.Reader(os.path.join(dir_path, 'external', 'ne_50m_coastline.shp'))
    x, y = [], []
    for shape in sf.shapes():
        px, py = zip(*shape.points)
        x += list(px) + [np.nan]
        y += list(py) + [np.nan]
    return np.array(x), np.array(y)

def BCADyear(year):
    '''Year as text in BC/AD.'''
    year = int(round(year))
    return f'{-year} BC' if year < 0 else f'{year} AD'      

def available_psvc():
    '''Print the available curves and models, with the elements and
    the time range read from their files.'''
    elements_by_columns = {7: 'D,I,F', 5: 'D,I', 3: 'F'}
    sections = [('local', 'Local PSVCs'), ('global', 'Global Models'),
                ('regional', 'Regional Models')]
    for kind, title in sections:
        print(f'\n{title}:')
        for name, c in catalog.items():
            if c['type'] != kind:
                continue
            if kind == 'local':
                m = np.atleast_2d(np.loadtxt(os.path.join(dir_path, 'curves', 'local', c['file'])))
                t, elements = m[:, 0], elements_by_columns.get(m.shape[1], '?')
            elif kind == 'global':
                t = np.loadtxt(os.path.join(dir_path, 'curves', 'gmodel', c['files'][2]))
                elements = 'D,I,F'
            else:
                t = np.loadtxt(os.path.join(dir_path, 'curves', 'rmodel', c['file']),
                               usecols=0, skiprows=1)
                elements = 'D,I,F'
            line = f'- {name}'
            if 'reference' in c:
                line += f" ({c['reference']})"
            line += f': [{elements}] [{BCADyear(t.min())} - {BCADyear(t.max())}]'
            if 'note' in c:
                line += f" ({c['note']})"
            print(line)

def log_sinh(x):
    '''log(sinh(x)) for x>0, without overflow for latge x'''
    return x + np.log1p(-np.exp(-2*x)) - np.log(2)

def cvp_direction(D, I, lat, lon, nlat, nlon):
    '''Relocate directions (D, I) from (lat, lon) to (nlat, nlon) through the
    virtual geomagnetic pole (CVP). D and I can be numbers or arrays.'''
    D, I = radians(D), radians(I)
    lat, lon, nlat, nlon = radians([lat, lon, nlat, nlon])
    clat, nclat = pi / 2 - lat, pi / 2 - nlat
    # dipolar coefficients
    g10 = - (cos(D) * cos(I) * sin(clat) + 0.5 * sin(I) * cos(clat))
    g11 = cos(D) * cos(I) * cos(clat) * cos(lon) + sin(D) * cos(I) * sin(lon) - 0.5 * sin(I) * sin(clat) * cos(lon)
    h11 = - (-cos(D) * cos(I) * cos(clat) * sin(lon) + sin(D) * cos(I) * cos(lon) + 0.5 * sin(I) * sin(clat) * sin(lon))
    # x, y, z at the new site
    x = -g10 * sin(nclat) + cos(nclat) * (g11 * cos(nlon) + h11 * sin(nlon))
    y = g11 * sin(nlon) - h11 * cos(nlon)
    z = -2 * (g10 * cos(nclat) + sin(nclat) * (g11 * cos(nlon) + h11 * sin(nlon)))
    return cart2dir(x, y, z)


def vdm_intensity(F, I, nI):
    '''Relocate intensities with the virtual dipole moment, from the site with
    inclination I to the site with inclination nI.'''
    clat = arctan(2 / tan(radians(I)))
    nclat = arctan(2 / tan(radians(nI)))
    return F * sqrt((1 + 3 * cos(nclat) ** 2) / (1 + 3 * cos(clat) ** 2))


def vadm_intensity(F, lat, nlat):
    '''Relocate intensities with the virtual axial dipole moment, from
    latitude lat to latitude nlat.'''
    clat, nclat = pi / 2 - radians(lat), pi / 2 - radians(nlat)
    return F * sqrt((1 + 3 * cos(nclat) ** 2) / (1 + 3 * cos(clat) ** 2))

def relocation_bias(model, lat, lon, nlat, nlon):
    '''Error of the dipolar relocation from (lat, lon) to (nlat, nlon),
    predicted by a global model at each of its epochs:

        bias = model at (nlat, nlon) - model at (lat, lon) relocated to (nlat, nlon)

    with CVP for the direction and VDM for the intensity.
    Returns the epochs of the model and the biases bD, bI, bF.'''
    m = Model(model)
    D0, I0, F0, *_ = m.psvc(lat, lon)
    D1, I1, F1, *_ = m.psvc(nlat, nlon)

    Dr, Ir = cvp_direction(D0, I0, lat, lon, nlat, nlon)
    Fr = vdm_intensity(F0, I0, Ir)

    bD = (D1 - Dr + 180) % 360 - 180   # angular difference
    bI = I1 - Ir
    bF = F1 - Fr
    return m.t, bD, bI, bF