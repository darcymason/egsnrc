# From https://github.com/numba/numba/issues/4037#issuecomment-1365116652
import numba as nb
from dataclasses import dataclass

def jitdataclass(cls=None, *, extra_spec=[]):
    """
    Helper decorator to make it easier to numba jitclass dataclasses

    Inspired by https://github.com/numba/numba/issues/4037#issuecomment-907523015
    """
    def _jitdataclass(cls):
        dc_cls = dataclass(cls, eq=False, match_args=False) # match_args Python >= 3.10
        del dc_cls.__dataclass_params__
        del dc_cls.__dataclass_fields__
        return nb.experimental.jitclass(dc_cls, spec=extra_spec)

    if cls is not None:
        # We've been called without additional args - invoke actual decorator
        return _jitdataclass(cls)
    # We've been called with additional args - so return actual decorator
    # which Python calls for us
    return _jitdataclass

@jitdataclass
class Particle:
    status: nb.int32
    region: nb.int32
    energy: nb.float32
    z: nb.float32