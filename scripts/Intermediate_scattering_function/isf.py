import numpy as np
from .base import AnalysisBase
from ..core import groups


class ISF(AnalysisBase):
 
    def __init__(self, u, select='all', **kwargs):
        
        if isinstance(u, groups.UpdatingAtomGroup):
            raise TypeError("UpdatingAtomGroups are not valid for MSD "
                            "computation")

        super(ISF, self).__init__(u.universe.trajectory, **kwargs)

        # args
        self.select = select
    
        # local
        self.ag = u.select_atoms(self.select)
        self.n_particles = len(self.ag)
        self._position_array = None

        # result
        self.results.timeseries = None

    def _prepare(self):
        # self.n_frames only available here
        # these need to be zeroed prior to each run() call
        self.results.timeseries = np.zeros(self.n_frames)
        self._position_array = np.zeros(
            (self.n_frames, self.n_particles, 3))
        # self.results.timeseries not set here


    def _single_frame(self):
        r""" Constructs array of positions for MSD calculation.

        """
        # shape of position array set here, use span in last dimension
        # from this point on
        keys ={'xyz': [0, 1, 2]}
        self._dim = keys['xyz']
        self._position_array[self._frame_index] = (
            self.ag.positions[:, self._dim])

    def _conclude(self):
        r""" Calculates the MSD via the simple "windowed" algorithm.

        """
        lagtimes = np.arange(1, self.n_frames)
        k_vec = 2*np.array([1,1,1])
        positions = self._position_array.astype(np.float64)
        for lag in lagtimes:
            disp = positions[lag:, :, :] - positions[:-lag, :, :] #r(t)-r(t0)
            dr_mul_k = np.sum(disp*k_vec, axis=2)
            sinA_div_A = np.sin(dr_mul_k)/dr_mul_k
            self.results.timeseries[lag] = np.mean(sinA_div_A)

