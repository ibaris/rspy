# # -*- coding: utf-8 -*-
# """
# This code is from Gomez Dans: https://github.com/jgomezdans/brdf_kernels
# """
#
# from sys import exit
#
# import numpy
#
#
# class Kernels:
#     '''
#     Linear kernel models
#     '''
#
#     def __init__(self, vza, sza, raa, critical=1, RossHS=True, RecipFlag=True, HB=2.0, BR=1.0, MODISSPARSE=True,
#                  MODISDENSE=False, RossType='Thick', normalise=1, normalize=0, LiType='Transit', doIntegrals=True,
#                  BSAangles=[], nbar=45.0):
#         '''
#         The class creator sets up the kernels for some angle set. Default Li is MODISSPARSE parameter set
# 	The kernels are accessible from:
# 		self.Isotropic
# 		self.Ross
# 		self.Li
# 	The angles are accesible from:
# 		self.vza (or self.vzaDegrees)
# 		self.sza (or self.szaDegrees)
# 		self.raa (or self.raaDegrees)
# 		N.B. Hot spot direction is vza == sza and raa = 0.0
# 	Kernels integrals are acessible from:
# 		self.BSAangles (angles in degrees)
# 		self.BSA_Isotropic (directional-hemispherical integral of self.Isotropic)
# 		self.BSA_Ross (directional-hemispherical integral of self.Ross)
# 		self.BSA_Li (directional-hemispherical integral of self.Li)
# 		self.WSA_Isotropic (bi-hemispherical integral of self.Isotropic)
# 		self.WSA_Ross (bi-hemispherical integral of self.Ross)
# 		self.WSA_Li (bi-hemispherical integral of self.Li)
# 		N.B. You need to set the doIntegrals flag to True on creating an instance of the kernels class if you
# 		want access to integrals. The processing takes a bit of time.
# 	Printing methods are available:
# 		self.printIntegrals(header=True,reflectance=False)
# 		self.printKernels(header=True,reflectance=False)
#
# 	Required parameters:
#
#         @param vza: an array containg view zenith angles in degrees
#         @param sza: an array containing solar zenith angles in degrees
#         @param raa: an array containg relative azimuth angles in degrees
#
# 	Options:
#         @option critical=1: set to 1 to exit on error, 0 not to
#         @option RecipFlag=True: Li reciprocal flag
#         @option HB: Li kernel parameter HB
#         @option BR: Li kernel parameter
#         @option MODISSPARSE: set to True for default MODIS Li Sparse parameters (overrides BR and HB to 2.0 and 1.0)
#         @option MODISDENSE: set to True for default MODIS Li Dense parameters (override BR and HB to 2.0 and 2.5)
#         @option RossType: set to 'Thin' for Ross Thin (default) else 'Thick'
#         @option LiType: set to 'Sparse' for LiSparse (default). Other options: 'Roujean', 'Dense'
#         @option normalise: set to 1 to make kernels 0 at nadir view illumination (default), set to 0 for no normalisation (can also use US spelling, i.e. normalize)
#         @option doIntegrals: set to True to calculate integrals of kernels numerically. Set to False not to calculate them. At some point will have Approx flag here as well.
#         @option BSAangles: solar zenith angles at which to calculate directional-hemispherical integral of kernels (default 0-89 in steps of 1 degree). Units: degrees.
#         @option nbar: the sza at which the isotropic term is set to if normalise=1 is turned on (default 0)
#
# 	Notes:
# 	Requires numpy. If you do integrals, this also requires scipy (or rather scipy.integrate)
# 	If you want to mimic the results in Wanner et al. 1995, I've set a special function called self.mimic at the end here.
#
#
#         '''
#         self.setup(critical=critical, RecipFlag=RecipFlag, RossHS=RossHS, HB=HB, BR=BR, MODISSPARSE=MODISSPARSE,
#                    MODISDENSE=MODISDENSE, RossType=RossType, normalise=normalise, normalize=normalize, LiType=LiType,
#                    doIntegrals=doIntegrals, BSAangles=BSAangles, nbar=nbar)
#         self.setAngleInfo(vza, sza, raa)
#         self.doKernels()
#         self.postProcess()
#
#     def setup(self, critical=1, RecipFlag=True, RossHS=True, HB=2.0, BR=1.0, MODISSPARSE=True, MODISDENSE=False,
#               RossType='Thick', normalise=1, normalize=0, LiType='Sparse', doIntegrals=True, BSAangles=[], nbar=0.0):
#         self.nbar = nbar
#         self.NEARLYZERO = 1e-20
#         self.critical = critical
#         self.FILE = -1
#         self.outputFile = ''
#         # kernel options etc.
#         self.LiType = LiType
#         self.RossHS = RossHS
#         self.doIntegrals = doIntegrals
#         if MODISDENSE == True:
#             LiType = 'Dense'
#             self.HB = 2.0
#             self.BR = 2.5
#         else:
#             if MODISSPARSE == True:
#                 LiType = 'Sparse'
#                 self.HB = 2.0
#                 self.BR = 1.0
#             else:
#                 self.HB = HB
#                 self.BR = BR
#         # self.LiType = LiType
#         self.RossType = RossType
#         self.normalise = normalise
#         self.RecipFlag = RecipFlag
#         # some useful numbers
#         self.M_PI = numpy.pi
#         self.M_PI_2 = self.M_PI * 0.5
#         self.M_PI_4 = self.M_PI * 0.25
#         self.M_1_PI = 1.0 / self.M_PI
#
#         self.normalise = 0
#         self.integrateKernels(BSAangles=BSAangles)
#
#         if (normalise >= 1 or normalize >= 1):
#             self.normalise = max(normalise, normalize)
#
#     def postProcess(self):
#         '''
#         Private method for dealing with normalisation
#         '''
#         self.LiNorm = 0.
#         self.RossNorm = 0.
#         self.IsotropicNorm = 0.
#         # if we are normalising the last element of self.Isotropic, self.Ross and self.Li  contain the nadir-nadir kernel
#         if self.normalise >= 1:
#             # normalise nbar-nadir (so kernel is 0 at nbar-nadir)
#             self.RossNorm = self.Ross[-1]
#             self.LiNorm = self.Li[-1]
#             self.Ross = self.Ross - self.RossNorm
#             self.Li = self.Li - self.LiNorm
#             # depreciate length of arrays (well, teh ones we'll use again in any case)
#             self.Ross = self.Ross[0:-1]
#             self.Li = self.Li[0:-1]
#             self.Isotropic = self.Isotropic[0:-1]
#             self.vzaDegrees = self.vzaDegrees[0:-1]
#             self.szaDegrees = self.szaDegrees[0:-1]
#             self.raaDegrees = self.raaDegrees[0:-1]
#             self.N = len(self.vzaDegrees)
#             self.vza = self.vza[0:-1]
#             self.sza = self.sza[0:-1]
#             self.raa = self.raa[0:-1]
#
#     def doKernels(self):
#         '''
#         Private method to run the various kernel methods
#         '''
#         # the kernels
#         self.IsotropicKernel()
#         self.RossKernel()
#         self.LiKernel()
#
#     def setAngleInfo(self, vza, sza, raa):
#         '''
#         Private method to store and organise the input angle data
#         '''
#         self.vzaDegrees = numpy.array([vza]).flatten()
#         self.szaDegrees = numpy.array([sza]).flatten()
#         self.raaDegrees = numpy.array([raa]).flatten()
#         self.N = len(self.vzaDegrees)
#
#         if (self.N != len(self.szaDegrees) or self.N != len(self.raaDegrees)):
#             self.error('kernels: inconsistent number of samples in vza, sza and raa data: ' + str(
#                 len(self.vzaDegrees)) + ', ' + str(len(self.szaDegrees)) + ', ' + str(len(self.raaDegrees)),
#                        critical=self.critical)
#             print self.vzaDegrees
#             print self.szaDegrees
#             print self.raaDegrees
#             return [-1]
#
#         if (self.normalise >= 1):
#             # calculate nadir term by extending array
#             self.vzaDegrees = numpy.array(list(self.vzaDegrees) + [0.0]).flatten()
#             self.szaDegrees = numpy.array(list(self.szaDegrees) + [self.nbar]).flatten()
#             self.raaDegrees = numpy.array(list(self.raaDegrees) + [0.0]).flatten()
#             # not N is one too many now
#             self.N = len(self.vzaDegrees)
#
#         self.vza = self.dtor(self.vzaDegrees)
#         self.sza = self.dtor(self.szaDegrees)  # -1 to make HS direction for raa = 0
#         self.raa = self.dtor(self.raaDegrees)
#         w = numpy.where(self.vza < 0)[0]
#         self.vza[w] = -self.vza[w]
#         self.raa[w] = self.raa[w] + numpy.pi
#         w = numpy.where(self.sza < 0)[0]
#         self.sza[w] = -self.sza[w]
#         self.raa[w] = self.raa[w] + numpy.pi
#
#     def integrateKernels(self, BSAangles=[]):
#         '''
#         Private method to call integration functions for the kernels
#
#
#              NB - this overwrites all kernel info ... so be careful how/where you call it
#             @option: BSAangles=[] allows the user to set the sza angles at which directional-hemispherical intergal is calculated, else steps of 1 degree from 0 to 89 (though I wouldnt trust it down to 90)
#             This function can be rather slow, so using fewer samples or an approximate function may be a god idea
#         '''
#         if (self.doIntegrals == False):
#             return;
#
#         import scipy.integrate
#         if BSAangles == []:
#             BSAangles = numpy.array(range(90)) * 1.0
#
#         self.BSAangles = numpy.array(BSAangles).flatten()
#
#         # isotropic integral
#         self.BSA_Isotropic = numpy.zeros(len(self.BSAangles)) + 1.0
#         self.BSA_Ross = numpy.zeros(len(self.BSAangles))
#         self.BSA_Li = numpy.zeros(len(self.BSAangles))
#         self.BSA_Isotropic_error = numpy.zeros(len(self.BSAangles))
#         self.BSA_Ross_error = numpy.zeros(len(self.BSAangles))
#         self.BSA_Li_error = numpy.zeros(len(self.BSAangles))
#
#         i = 0
#         mu = numpy.cos(self.BSAangles * numpy.pi / 180.)
#         for sza in self.BSAangles:
#             # ross integral
#             self.BSA_Ross[i], self.BSA_Ross_error[i] = scipy.integrate.dblquad(RossFunctionForIntegral, 0.0, 1.0,
#                                                                                gfun, hfun, args=(sza, self))
#             self.BSA_Li[i], self.BSA_Li_error[i] = scipy.integrate.dblquad(LiFunctionForIntegral, 0.0, 1.0, gfun,
#                                                                            hfun, args=(sza, self))
#             i = i + 1
#         self.WSA_Ross = -2.0 * scipy.integrate.simps(self.BSA_Ross * mu, mu)
#         self.WSA_Li = -2.0 * scipy.integrate.simps(self.BSA_Li * mu, mu)
#         return
#
#     def GetPhaang(self):
#         '''
#         Private method to calculate Phase angle component of kernel
#         '''
#         self.cosphaang = self.cos1 * self.cos2 + self.sin1 * self.sin2 * self.cos3
#         # better check the bounds before arccos ... just to be safe
#         w = numpy.where(self.cosphaang < -1)[0]
#         self.cosphaang[w] = -1.0
#         w = numpy.where(self.cosphaang > 1)[0]
#         self.cosphaang[w] = 1.0
#         self.phaang = numpy.arccos(self.cosphaang)
#         self.sinphaang = numpy.sin(self.phaang)
#         return
#
#     def RossKernelPart(self):
#         '''
# 	Private method to calculate main part of Ross kernel
# 	'''
#         self.cos1 = numpy.cos(self.vza)
#         self.cos2 = numpy.cos(self.sza)
#
#         self.sin1 = numpy.sin(self.vza)
#         self.sin2 = numpy.sin(self.sza)
#         self.cos3 = numpy.cos(self.raa)
#         self.GetPhaang()
#         self.rosselement = (self.M_PI_2 - self.phaang) * self.cosphaang + self.sinphaang
#         return
#
#     def GetDistance(self):
#         '''
#         Private method to get distance component of Li kernels
#         '''
#         temp = self.tan1 * self.tan1 + self.tan2 * self.tan2 - 2. * self.tan1 * self.tan2 * self.cos3;
#         w = numpy.where(temp < 0)[0]
#         temp[w] = 0.0
#         self.temp = temp  # used by other functions ??
#         distance = numpy.sqrt(temp)
#         return distance
#
#     def GetpAngles(self, tan1):
#         '''
#         Private method to do B/R transformation for ellipse shape
#         '''
#         t = self.BR * tan1
#         w = numpy.where(t < 0.)[0]
#         t[w] = 0.0
#         angp = numpy.arctan(t)
#         s = numpy.sin(angp)
#         c = numpy.cos(angp)
#         # have to make sure c isnt 0
#         w = numpy.where(c == 0)[0]
#         c[w] = self.NEARLYZERO
#         return c, s, t
#
#     def GetOverlap(self):
#         '''
#         Private method to do HB ratio transformation
#         '''
#         self.temp = 1. / self.cos1 + 1. / self.cos2
#
#         self.cost = self.HB * numpy.sqrt(
#             self.distance * self.distance + self.tan1 * self.tan1 * self.tan2 * self.tan2 * self.sin3 * self.sin3) / self.temp;
#         w = numpy.where(self.cost < -1)[0]
#         self.cost[w] = -1.0
#         w = numpy.where(self.cost > 1.0)[0]
#         self.cost[w] = 1.0
#         self.tvar = numpy.arccos(self.cost)
#         self.sint = numpy.sin(self.tvar)
#         self.overlap = self.M_1_PI * (self.tvar - self.sint * self.cost) * self.temp
#         w = numpy.where(self.overlap < 0)[0]
#         self.overlap[w] = 0.0
#         return
#
#     def RoujeanKernel(self):
#         '''
#         Private method - call to calculate Roujean shadowing kernel
#         '''
#         # first make sure its in range 0 to 2 pi
#         self.phi = numpy.abs((self.raa % (2. * numpy.pi)))
#         self.cos3 = numpy.cos(self.phi)
#         self.sin3 = numpy.sin(self.phi)
#         self.tan1 = numpy.tan(self.sza)
#         self.tan2 = numpy.tan(self.vza)
#
#         self.distance = self.GetDistance()
#         self.Li = 0.5 * self.M_1_PI * ((
#                                                self.M_PI - self.phi) * self.cos3 + self.sin3) * self.tan1 * self.tan2 - self.M_1_PI * (
#                           self.tan1 + self.tan2 + self.distance);
#         return
#
#     def LiKernel(self):
#         '''
#         Private method - call to calculate Li Kernel
#         '''
#         # at some point add in LiGround kernel & LiTransit
#         if self.LiType == 'Roujean':
#             return self.RoujeanKernel()
#         # first make sure its in range 0 to 2 pi
#         self.phi = numpy.abs((self.raa % (2. * numpy.pi)))
#         self.cos3 = numpy.cos(self.phi)
#         self.sin3 = numpy.sin(self.phi)
#         self.tanti = numpy.tan(self.sza)
#         self.tantv = numpy.tan(self.vza)
#         self.cos1, self.sin1, self.tan1 = self.GetpAngles(self.tantv);
#         self.cos2, self.sin2, self.tan2 = self.GetpAngles(self.tanti);
#         self.GetPhaang();  # sets cos & sin phase angle terms
#         self.distance = self.GetDistance();  # sets self.temp
#         self.GetOverlap();  # also sets self.temp
#         if self.LiType == 'Sparse':
#             if self.RecipFlag == True:
#                 self.Li = self.overlap - self.temp + 0.5 * (1. + self.cosphaang) / self.cos1 / self.cos2;
#             else:
#                 self.Li = self.overlap - self.temp + 0.5 * (1. + self.cosphaang) / self.cos1;
#         else:
#             if self.LiType == 'Dense':
#                 if self.RecipFlag:
#                     self.Li = (1.0 + self.cosphaang) / (
#                             self.cos1 * self.cos2 * (self.temp - self.overlap)) - 2.0;
#                 else:
#                     self.Li = (1.0 + self.cosphaang) / (self.cos1 * (self.temp - self.overlap)) - 2.0;
#             else:
#                 B = self.temp - self.overlap
#                 w = numpy.where(B <= 2.0)
#                 self.Li = B * 0.0
#                 if self.RecipFlag == True:
#                     Li = self.overlap - self.temp + 0.5 * (1. + self.cosphaang) / self.cos1 / self.cos2;
#                 else:
#                     Li = self.overlap - self.temp + 0.5 * (1. + self.cosphaang) / self.cos1;
#                 self.Li[w] = Li[w]
#
#                 w = numpy.where(B > 2.0)
#                 if self.RecipFlag:
#                     Li = (1.0 + self.cosphaang) / (self.cos1 * self.cos2 * (self.temp - self.overlap)) - 2.0;
#                 else:
#                     Li = (1.0 + self.cosphaang) / (self.cos1 * (self.temp - self.overlap)) - 2.0;
#                 self.Li[w] = Li[w]
#         return
#
#     def IsotropicKernel(self):
#         '''
#         Public method - call to calculate Isotropic kernel
#         '''
#         # default behaviour
#         self.Isotropic = numpy.zeros(self.N) + 1.0
#         return
#
#     def RossThin(self):
#         '''
#         Public method - call to calculate RossThin kernel
#         '''
#         self.RossKernelPart()
#         self.rosselement = self.rosselement / (self.cos1 * self.cos2)
#         return;
#
#     def RossThick(self):
#         '''
#         Public method - call to calculate RossThick kernel
#         '''
#         self.RossKernelPart()
#         self.rosselement = self.rosselement / (self.cos1 + self.cos2)
#         return;
#
#     def RossKernel(self):
#         '''
#         Public method - call to calculate Ross Kernel
#         '''
#         if self.RossType == 'Thin':
#             self.RossThin()
#         else:
#             self.RossThick()
#         self.Ross = self.rosselement
#         if self.RossHS != False:
#             if self.RossHS == True:
#                 self.RossHS = 0.25
#             self.Ross = self.Ross * (1 + 1 / (1 + self.phaang / self.RossHS))
#
#     def dtor(self, x):
#         '''
#         Public method to convert degrees to radians
#         '''
#         return x * numpy.pi / 180.0
#
#     def rtod(self, x):
#         '''
#         Public method to convert radians to degrees
#         '''
#         return x * 180. / numpy.pi
#
#     def error(self, msg, critical=0, newline=1, code=-1):
#         '''
#         Public method to do Class error reporting
#         @param msg: error message
#         @param critical: set to 1 if require exit (default critical=0)
#         @param newline: set to 0 if newline not required (default newline=0)
#         @param code: error code reported on exit if critical error (default code=-1)
#         '''
#         if newline == 1:
#             nl = '\n'
#         else:
#             nl = ''
#         print msg + nl
#         if critical == 1:
#             exit([code])
#
#     def printIntegrals(self, header=True, reflectance=False):
#         '''
#         Public method to print kernel integrals (to stdout only at present)
#         '''
#         if (header == True):
#             self.printer(
#                 '# ' + str(self.N) + ' samples Ross: ' + self.RossType + ' Li: ' + self.LiType + ' Reciprocal: ' + str(
#                     self.RecipFlag) + ' normalisation: ' + str(self.normalise) + ' HB ' + str(self.HB) + ' BR ' + str(
#                     self.BR) + '\n');
#             self.printer('# WSA: Isotropic 1.0 Ross ' + str(self.WSA_Ross) + ' Li ' + str(self.WSA_Li))
#             self.printer('# 1: SZA (degrees) 2: BSA Isotropic 3: BSA Ross 4: BSA Li')
#             if (reflectance == True):
#                 self.printer(' ');
#             self.printer('\n');
#
#         for i in range(len(self.BSAangles)):
#             self.printer(
#                 str(self.BSAangles[i]) + ' ' + str(self.BSA_Isotropic[i]) + ' ' + str(self.BSA_Ross[i]) + ' ' + str(
#                     self.BSA_Li[i]))
#             # print refl data if wanted
#             if (reflectance == True):
#                 self.printer(' ');
#             self.printer('\n');
#         return
#
#     def printKernels(self, header=True, reflectance=False, file=False):
#         '''
#         Public method to print kernel values (to stdout only at present)
# 	'''
#         if (file != False):
#             if (file != self.outputFile and self.FILE != -1):
#                 self.FILE.close()
#             self.outputFile = file
#             self.FILE = open(self.outputFile, 'w')
#
#         if (header == True):
#             self.printer(
#                 '# ' + str(self.N) + ' samples Ross: ' + self.RossType + ' Li: ' + self.LiType + ' Reciprocal: ' + str(
#                     self.RecipFlag) + ' normalisation: ' + str(self.normalise) + ' HB ' + str(self.HB) + ' BR ' + str(
#                     self.BR) + '\n');
#             self.printer('# 1: VZA (degrees) 2: SZA (degrees) 3: RAA (degrees) 4: Isotropic 5: Ross 6: Li')
#             if (reflectance == True):
#                 self.printer(' ');
#             self.printer('\n');
#
#         for i in range(self.N):
#             self.printer(
#                 str(self.vzaDegrees[i]) + ' ' + str(self.szaDegrees[i]) + ' ' + str(self.raaDegrees[i]) + ' ' + str(
#                     self.Isotropic[i]) + ' ' + str(self.Ross[i]) + ' ' + str(self.Li[i]))
#             # print refl data if wanted
#             if (reflectance == True):
#                 self.printer(' ');
#             self.printer('\n');
#         return
#
#     def printer(self, msg):
#         '''
#         Public print method ... make more flexible eg for printing to files at some point
#         '''
#         if (self.FILE == -1):
#             print msg,
#         else:
#             self.FILE.write(msg)
#
#
# # some things required for the numerical integration
#
# def _Kernelsgfun(x):
#     return 0.0
#
#
# def _Kernelshfun(x):
#     return 2.0 * numpy.pi
#
#
# def RossFunctionForIntegral(phi, mu, sza, self):
#     # print phi
#     # print mu
#     # print sza
#     # print '========'
#     vza = numpy.arccos(mu)
#     raa = self.rtod(phi)
#     self.setAngleInfo(vza, sza, raa)
#     self.RossKernel()
#     return mu * self.Ross[0] / numpy.pi
#
#
# def LiFunctionForIntegral(phi, mu, sza, self):
#     # print phi
#     # print mu
#     # print sza
#     # print '========'
#     vza = numpy.arccos(mu)
#     raa = self.rtod(phi)
#     self.setAngleInfo(vza, sza, raa)
#     self.LiKernel()
#     return mu * self.Li[0] / numpy.pi
#
#
# import numpy as np
#
# n = 100
#
# iza, vza, iaa, vaa, nbar = np.random.uniform(1, np.pi / 2, n), np.random.uniform(1, np.pi / 2, n), \
#                            np.random.uniform(1, 2 * np.pi, n), np.random.uniform(1, 2 * np.pi, n), \
#                            np.random.uniform(1, np.pi / 4, 1)
#
# izaDeg, vzaDeg, iaaDeg, vaaDeg, nbarDeg = np.rad2deg(iza), np.rad2deg(vza), np.rad2deg(
#     iaa), np.rad2deg(vaa), np.rad2deg(nbar)
#
# raa = np.abs(iaa - vaa)
# raaDeg = iaaDeg - vaaDeg
#
# angles = np.array([iza, vza, raa])
# np.savetxt("/home/ibaris/Dropbox/GitHub/rspy/rspy/test/test_phase/data/angles.out", angles.transpose())
#
# k = Kernels(vzaDeg, izaDeg, raaDeg, doIntegrals=False, HB=2.0, BR=1.0)
# phaang = k.phaang[0:-1]
# sinphaang = k.sinphaang[0:-1]
# cosphaang = k.cosphaang[0:-1]
# O = k.overlap[0:-1]
# D = k.distance[0:-1]
# piza = k.tan2[0:-1]
# pvza = k.tan1[0:-1]
#
# phaangs = np.array([phaang, sinphaang, cosphaang, O, D, piza, pvza])
# np.savetxt("/home/ibaris/Dropbox/GitHub/rspy/rspy/test/test_phase/data/results_dans.out", phaangs.transpose())
