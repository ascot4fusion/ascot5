# -*- coding: utf-8 -*-
#
# TARGET arch is: ['-I/usr/include/hdf5/serial', '-I/usr/lib/gcc/x86_64-linux-gnu/11/include/']
# WORD_SIZE is: 8
# POINTER_SIZE is: 8
# LONGDOUBLE_SIZE is: 16
#
import ctypes


class AsDictMixin:
    @classmethod
    def as_dict(cls, self):
        result = {}
        if not isinstance(self, AsDictMixin):
            # not a structure, assume it's already a python object
            return self
        if not hasattr(cls, "_fields_"):
            return result
        # sys.version_info >= (3, 5)
        # for (field, *_) in cls._fields_:  # noqa
        for field_tuple in cls._fields_:  # noqa
            field = field_tuple[0]
            if field.startswith('PADDING_'):
                continue
            value = getattr(self, field)
            type_ = type(value)
            if hasattr(value, "_length_") and hasattr(value, "_type_"):
                # array
                if not hasattr(type_, "as_dict"):
                    value = [v for v in value]
                else:
                    type_ = type_._type_
                    value = [type_.as_dict(v) for v in value]
            elif hasattr(value, "contents") and hasattr(value, "_type_"):
                # pointer
                try:
                    if not hasattr(type_, "as_dict"):
                        value = value.contents
                    else:
                        type_ = type_._type_
                        value = type_.as_dict(value.contents)
                except ValueError:
                    # nullptr
                    value = None
            elif isinstance(value, AsDictMixin):
                # other structure
                value = type_.as_dict(value)
            result[field] = value
        return result


class Structure(ctypes.Structure, AsDictMixin):

    def __init__(self, *args, **kwds):
        # We don't want to use positional arguments fill PADDING_* fields

        args = dict(zip(self.__class__._field_names_(), args))
        args.update(kwds)
        super(Structure, self).__init__(**args)

    @classmethod
    def _field_names_(cls):
        if hasattr(cls, '_fields_'):
            return (f[0] for f in cls._fields_ if not f[0].startswith('PADDING'))
        else:
            return ()

    @classmethod
    def get_type(cls, field):
        for f in cls._fields_:
            if f[0] == field:
                return f[1]
        return None

    @classmethod
    def bind(cls, bound_fields):
        fields = {}
        for name, type_ in cls._fields_:
            if hasattr(type_, "restype"):
                if name in bound_fields:
                    if bound_fields[name] is None:
                        fields[name] = type_()
                    else:
                        # use a closure to capture the callback from the loop scope
                        fields[name] = (
                            type_((lambda callback: lambda *args: callback(*args))(
                                bound_fields[name]))
                        )
                    del bound_fields[name]
                else:
                    # default callback implementation (does nothing)
                    try:
                        default_ = type_(0).restype().value
                    except TypeError:
                        default_ = None
                    fields[name] = type_((
                        lambda default_: lambda *args: default_)(default_))
            else:
                # not a callback function, use default initialization
                if name in bound_fields:
                    fields[name] = bound_fields[name]
                    del bound_fields[name]
                else:
                    fields[name] = type_()
        if len(bound_fields) != 0:
            raise ValueError(
                "Cannot bind the following unknown callback(s) {}.{}".format(
                    cls.__name__, bound_fields.keys()
            ))
        return cls(**fields)


class Union(ctypes.Union, AsDictMixin):
    pass



_libraries = {}
# Ugly solution to find libascot.so
from pathlib import Path
libpath = str(Path(__file__).absolute().parent.parent.parent) \
+ "/build/libascot.so"
_libraries['libascot.so'] = ctypes.CDLL(libpath)
c_int128 = ctypes.c_ubyte*16
c_uint128 = c_int128
void = None
if ctypes.sizeof(ctypes.c_longdouble) == 16:
    c_long_double_t = ctypes.c_longdouble
else:
    c_long_double_t = ctypes.c_ubyte*16

def string_cast(char_pointer, encoding='utf-8', errors='strict'):
    value = ctypes.cast(char_pointer, ctypes.c_char_p).value
    if value is not None and encoding is not None:
        value = value.decode(encoding, errors=errors)
    return value


def char_pointer_cast(string, encoding='utf-8'):
    if encoding is not None:
        try:
            string = string.encode(encoding)
        except AttributeError:
            # In Python3, bytes has no encode attribute
            pass
    string = ctypes.c_char_p(string)
    return ctypes.cast(string, ctypes.POINTER(ctypes.c_char))





integer = ctypes.c_int64
real = ctypes.c_double
class struct_c__SA_B_GS_offload_data(Structure):
    pass

struct_c__SA_B_GS_offload_data._pack_ = 1 # source:False
struct_c__SA_B_GS_offload_data._fields_ = [
    ('R0', ctypes.c_double),
    ('z0', ctypes.c_double),
    ('raxis', ctypes.c_double),
    ('zaxis', ctypes.c_double),
    ('B_phi0', ctypes.c_double),
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('psi_mult', ctypes.c_double),
    ('psi_coeff', ctypes.c_double * 13),
    ('Nripple', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('a0', ctypes.c_double),
    ('alpha0', ctypes.c_double),
    ('delta0', ctypes.c_double),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

B_GS_offload_data = struct_c__SA_B_GS_offload_data
class struct_c__SA_B_GS_data(Structure):
    pass

struct_c__SA_B_GS_data._pack_ = 1 # source:False
struct_c__SA_B_GS_data._fields_ = [
    ('R0', ctypes.c_double),
    ('z0', ctypes.c_double),
    ('raxis', ctypes.c_double),
    ('zaxis', ctypes.c_double),
    ('B_phi0', ctypes.c_double),
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('psi_mult', ctypes.c_double),
    ('psi_coeff', ctypes.POINTER(ctypes.c_double)),
    ('Nripple', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('a0', ctypes.c_double),
    ('alpha0', ctypes.c_double),
    ('delta0', ctypes.c_double),
]

B_GS_data = struct_c__SA_B_GS_data
B_GS_init_offload = _libraries['libascot.so'].B_GS_init_offload
B_GS_init_offload.restype = ctypes.c_int32
B_GS_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_GS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_GS_free_offload = _libraries['libascot.so'].B_GS_free_offload
B_GS_free_offload.restype = None
B_GS_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_GS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_GS_init = _libraries['libascot.so'].B_GS_init
B_GS_init.restype = None
B_GS_init.argtypes = [ctypes.POINTER(struct_c__SA_B_GS_data), ctypes.POINTER(struct_c__SA_B_GS_offload_data), ctypes.POINTER(ctypes.c_double)]
a5err = ctypes.c_uint64
B_GS_eval_B = _libraries['libascot.so'].B_GS_eval_B
B_GS_eval_B.restype = a5err
B_GS_eval_B.argtypes = [ctypes.c_double * 3, real, real, real, ctypes.POINTER(struct_c__SA_B_GS_data)]
B_GS_eval_psi = _libraries['libascot.so'].B_GS_eval_psi
B_GS_eval_psi.restype = a5err
B_GS_eval_psi.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, ctypes.POINTER(struct_c__SA_B_GS_data)]
B_GS_eval_psi_dpsi = _libraries['libascot.so'].B_GS_eval_psi_dpsi
B_GS_eval_psi_dpsi.restype = a5err
B_GS_eval_psi_dpsi.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_GS_data)]
B_GS_eval_rho_drho = _libraries['libascot.so'].B_GS_eval_rho_drho
B_GS_eval_rho_drho.restype = a5err
B_GS_eval_rho_drho.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_GS_data)]
B_GS_eval_B_dB = _libraries['libascot.so'].B_GS_eval_B_dB
B_GS_eval_B_dB.restype = a5err
B_GS_eval_B_dB.argtypes = [ctypes.c_double * 12, real, real, real, ctypes.POINTER(struct_c__SA_B_GS_data)]
B_GS_get_axis_rz = _libraries['libascot.so'].B_GS_get_axis_rz
B_GS_get_axis_rz.restype = a5err
B_GS_get_axis_rz.argtypes = [ctypes.c_double * 2, ctypes.POINTER(struct_c__SA_B_GS_data)]
class struct_c__SA_B_2DS_offload_data(Structure):
    pass

struct_c__SA_B_2DS_offload_data._pack_ = 1 # source:False
struct_c__SA_B_2DS_offload_data._fields_ = [
    ('n_r', ctypes.c_int32),
    ('n_z', ctypes.c_int32),
    ('r_min', ctypes.c_double),
    ('r_max', ctypes.c_double),
    ('z_min', ctypes.c_double),
    ('z_max', ctypes.c_double),
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('axis_r', ctypes.c_double),
    ('axis_z', ctypes.c_double),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
]

B_2DS_offload_data = struct_c__SA_B_2DS_offload_data
class struct_c__SA_B_2DS_data(Structure):
    pass

class struct_c__SA_interp2D_data(Structure):
    pass

struct_c__SA_interp2D_data._pack_ = 1 # source:False
struct_c__SA_interp2D_data._fields_ = [
    ('n_x', ctypes.c_int32),
    ('n_y', ctypes.c_int32),
    ('bc_x', ctypes.c_int32),
    ('bc_y', ctypes.c_int32),
    ('x_min', ctypes.c_double),
    ('x_max', ctypes.c_double),
    ('x_grid', ctypes.c_double),
    ('y_min', ctypes.c_double),
    ('y_max', ctypes.c_double),
    ('y_grid', ctypes.c_double),
    ('c', ctypes.POINTER(ctypes.c_double)),
]

struct_c__SA_B_2DS_data._pack_ = 1 # source:False
struct_c__SA_B_2DS_data._fields_ = [
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('axis_r', ctypes.c_double),
    ('axis_z', ctypes.c_double),
    ('psi', struct_c__SA_interp2D_data),
    ('B_r', struct_c__SA_interp2D_data),
    ('B_phi', struct_c__SA_interp2D_data),
    ('B_z', struct_c__SA_interp2D_data),
]

B_2DS_data = struct_c__SA_B_2DS_data
B_2DS_init_offload = _libraries['libascot.so'].B_2DS_init_offload
B_2DS_init_offload.restype = ctypes.c_int32
B_2DS_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_2DS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_2DS_free_offload = _libraries['libascot.so'].B_2DS_free_offload
B_2DS_free_offload.restype = None
B_2DS_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_2DS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_2DS_init = _libraries['libascot.so'].B_2DS_init
B_2DS_init.restype = None
B_2DS_init.argtypes = [ctypes.POINTER(struct_c__SA_B_2DS_data), ctypes.POINTER(struct_c__SA_B_2DS_offload_data), ctypes.POINTER(ctypes.c_double)]
B_2DS_eval_psi = _libraries['libascot.so'].B_2DS_eval_psi
B_2DS_eval_psi.restype = a5err
B_2DS_eval_psi.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, ctypes.POINTER(struct_c__SA_B_2DS_data)]
B_2DS_eval_psi_dpsi = _libraries['libascot.so'].B_2DS_eval_psi_dpsi
B_2DS_eval_psi_dpsi.restype = a5err
B_2DS_eval_psi_dpsi.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_2DS_data)]
B_2DS_eval_rho_drho = _libraries['libascot.so'].B_2DS_eval_rho_drho
B_2DS_eval_rho_drho.restype = a5err
B_2DS_eval_rho_drho.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_2DS_data)]
B_2DS_eval_B = _libraries['libascot.so'].B_2DS_eval_B
B_2DS_eval_B.restype = a5err
B_2DS_eval_B.argtypes = [ctypes.c_double * 3, real, real, real, ctypes.POINTER(struct_c__SA_B_2DS_data)]
B_2DS_eval_B_dB = _libraries['libascot.so'].B_2DS_eval_B_dB
B_2DS_eval_B_dB.restype = a5err
B_2DS_eval_B_dB.argtypes = [ctypes.c_double * 12, real, real, real, ctypes.POINTER(struct_c__SA_B_2DS_data)]
B_2DS_get_axis_rz = _libraries['libascot.so'].B_2DS_get_axis_rz
B_2DS_get_axis_rz.restype = a5err
B_2DS_get_axis_rz.argtypes = [ctypes.c_double * 2, ctypes.POINTER(struct_c__SA_B_2DS_data)]
class struct_c__SA_B_3DS_offload_data(Structure):
    pass

struct_c__SA_B_3DS_offload_data._pack_ = 1 # source:False
struct_c__SA_B_3DS_offload_data._fields_ = [
    ('psigrid_n_r', ctypes.c_int32),
    ('psigrid_n_z', ctypes.c_int32),
    ('psigrid_r_min', ctypes.c_double),
    ('psigrid_r_max', ctypes.c_double),
    ('psigrid_z_min', ctypes.c_double),
    ('psigrid_z_max', ctypes.c_double),
    ('Bgrid_n_r', ctypes.c_int32),
    ('Bgrid_n_z', ctypes.c_int32),
    ('Bgrid_r_min', ctypes.c_double),
    ('Bgrid_r_max', ctypes.c_double),
    ('Bgrid_z_min', ctypes.c_double),
    ('Bgrid_z_max', ctypes.c_double),
    ('Bgrid_n_phi', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('Bgrid_phi_min', ctypes.c_double),
    ('Bgrid_phi_max', ctypes.c_double),
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('axis_r', ctypes.c_double),
    ('axis_z', ctypes.c_double),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

B_3DS_offload_data = struct_c__SA_B_3DS_offload_data
class struct_c__SA_B_3DS_data(Structure):
    pass

class struct_c__SA_interp3D_data(Structure):
    pass

struct_c__SA_interp3D_data._pack_ = 1 # source:False
struct_c__SA_interp3D_data._fields_ = [
    ('n_x', ctypes.c_int32),
    ('n_y', ctypes.c_int32),
    ('n_z', ctypes.c_int32),
    ('bc_x', ctypes.c_int32),
    ('bc_y', ctypes.c_int32),
    ('bc_z', ctypes.c_int32),
    ('x_min', ctypes.c_double),
    ('x_max', ctypes.c_double),
    ('x_grid', ctypes.c_double),
    ('y_min', ctypes.c_double),
    ('y_max', ctypes.c_double),
    ('y_grid', ctypes.c_double),
    ('z_min', ctypes.c_double),
    ('z_max', ctypes.c_double),
    ('z_grid', ctypes.c_double),
    ('c', ctypes.POINTER(ctypes.c_double)),
]

struct_c__SA_B_3DS_data._pack_ = 1 # source:False
struct_c__SA_B_3DS_data._fields_ = [
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('axis_r', ctypes.c_double),
    ('axis_z', ctypes.c_double),
    ('psi', struct_c__SA_interp2D_data),
    ('B_r', struct_c__SA_interp3D_data),
    ('B_phi', struct_c__SA_interp3D_data),
    ('B_z', struct_c__SA_interp3D_data),
]

B_3DS_data = struct_c__SA_B_3DS_data
B_3DS_init_offload = _libraries['libascot.so'].B_3DS_init_offload
B_3DS_init_offload.restype = ctypes.c_int32
B_3DS_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_3DS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_3DS_free_offload = _libraries['libascot.so'].B_3DS_free_offload
B_3DS_free_offload.restype = None
B_3DS_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_3DS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_3DS_init = _libraries['libascot.so'].B_3DS_init
B_3DS_init.restype = None
B_3DS_init.argtypes = [ctypes.POINTER(struct_c__SA_B_3DS_data), ctypes.POINTER(struct_c__SA_B_3DS_offload_data), ctypes.POINTER(ctypes.c_double)]
B_3DS_eval_psi = _libraries['libascot.so'].B_3DS_eval_psi
B_3DS_eval_psi.restype = a5err
B_3DS_eval_psi.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, ctypes.POINTER(struct_c__SA_B_3DS_data)]
B_3DS_eval_psi_dpsi = _libraries['libascot.so'].B_3DS_eval_psi_dpsi
B_3DS_eval_psi_dpsi.restype = a5err
B_3DS_eval_psi_dpsi.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_3DS_data)]
B_3DS_eval_rho_drho = _libraries['libascot.so'].B_3DS_eval_rho_drho
B_3DS_eval_rho_drho.restype = a5err
B_3DS_eval_rho_drho.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_3DS_data)]
B_3DS_eval_B = _libraries['libascot.so'].B_3DS_eval_B
B_3DS_eval_B.restype = a5err
B_3DS_eval_B.argtypes = [ctypes.c_double * 3, real, real, real, ctypes.POINTER(struct_c__SA_B_3DS_data)]
B_3DS_eval_B_dB = _libraries['libascot.so'].B_3DS_eval_B_dB
B_3DS_eval_B_dB.restype = a5err
B_3DS_eval_B_dB.argtypes = [ctypes.c_double * 12, real, real, real, ctypes.POINTER(struct_c__SA_B_3DS_data)]
B_3DS_get_axis_rz = _libraries['libascot.so'].B_3DS_get_axis_rz
B_3DS_get_axis_rz.restype = a5err
B_3DS_get_axis_rz.argtypes = [ctypes.c_double * 2, ctypes.POINTER(struct_c__SA_B_3DS_data)]
class struct_c__SA_B_STS_offload_data(Structure):
    pass

struct_c__SA_B_STS_offload_data._pack_ = 1 # source:False
struct_c__SA_B_STS_offload_data._fields_ = [
    ('psigrid_n_r', ctypes.c_int32),
    ('psigrid_n_z', ctypes.c_int32),
    ('psigrid_n_phi', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('psigrid_r_min', ctypes.c_double),
    ('psigrid_r_max', ctypes.c_double),
    ('psigrid_z_min', ctypes.c_double),
    ('psigrid_z_max', ctypes.c_double),
    ('psigrid_phi_min', ctypes.c_double),
    ('psigrid_phi_max', ctypes.c_double),
    ('Bgrid_n_r', ctypes.c_int32),
    ('Bgrid_n_z', ctypes.c_int32),
    ('Bgrid_n_phi', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('Bgrid_r_min', ctypes.c_double),
    ('Bgrid_r_max', ctypes.c_double),
    ('Bgrid_z_min', ctypes.c_double),
    ('Bgrid_z_max', ctypes.c_double),
    ('Bgrid_phi_min', ctypes.c_double),
    ('Bgrid_phi_max', ctypes.c_double),
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('offload_array_length', ctypes.c_int32),
    ('n_axis', ctypes.c_int32),
    ('axis_min', ctypes.c_double),
    ('axis_max', ctypes.c_double),
    ('axis_grid', ctypes.c_double),
]

B_STS_offload_data = struct_c__SA_B_STS_offload_data
class struct_c__SA_B_STS_data(Structure):
    pass

class struct_c__SA_linint1D_data(Structure):
    pass

struct_c__SA_linint1D_data._pack_ = 1 # source:False
struct_c__SA_linint1D_data._fields_ = [
    ('n_x', ctypes.c_int32),
    ('bc_x', ctypes.c_int32),
    ('x_min', ctypes.c_double),
    ('x_max', ctypes.c_double),
    ('x_grid', ctypes.c_double),
    ('c', ctypes.POINTER(ctypes.c_double)),
]

struct_c__SA_B_STS_data._pack_ = 1 # source:False
struct_c__SA_B_STS_data._fields_ = [
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('axis_r', struct_c__SA_linint1D_data),
    ('axis_z', struct_c__SA_linint1D_data),
    ('psi', struct_c__SA_interp3D_data),
    ('B_r', struct_c__SA_interp3D_data),
    ('B_phi', struct_c__SA_interp3D_data),
    ('B_z', struct_c__SA_interp3D_data),
]

B_STS_data = struct_c__SA_B_STS_data
B_STS_init_offload = _libraries['libascot.so'].B_STS_init_offload
B_STS_init_offload.restype = ctypes.c_int32
B_STS_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_STS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_STS_free_offload = _libraries['libascot.so'].B_STS_free_offload
B_STS_free_offload.restype = None
B_STS_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_STS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_STS_init = _libraries['libascot.so'].B_STS_init
B_STS_init.restype = None
B_STS_init.argtypes = [ctypes.POINTER(struct_c__SA_B_STS_data), ctypes.POINTER(struct_c__SA_B_STS_offload_data), ctypes.POINTER(ctypes.c_double)]
B_STS_eval_psi = _libraries['libascot.so'].B_STS_eval_psi
B_STS_eval_psi.restype = a5err
B_STS_eval_psi.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, ctypes.POINTER(struct_c__SA_B_STS_data)]
B_STS_eval_psi_dpsi = _libraries['libascot.so'].B_STS_eval_psi_dpsi
B_STS_eval_psi_dpsi.restype = a5err
B_STS_eval_psi_dpsi.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_STS_data)]
B_STS_eval_rho_drho = _libraries['libascot.so'].B_STS_eval_rho_drho
B_STS_eval_rho_drho.restype = a5err
B_STS_eval_rho_drho.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_STS_data)]
B_STS_eval_B = _libraries['libascot.so'].B_STS_eval_B
B_STS_eval_B.restype = a5err
B_STS_eval_B.argtypes = [ctypes.c_double * 3, real, real, real, ctypes.POINTER(struct_c__SA_B_STS_data)]
B_STS_eval_B_dB = _libraries['libascot.so'].B_STS_eval_B_dB
B_STS_eval_B_dB.restype = a5err
B_STS_eval_B_dB.argtypes = [ctypes.c_double * 12, real, real, real, ctypes.POINTER(struct_c__SA_B_STS_data)]
B_STS_get_axis_rz = _libraries['libascot.so'].B_STS_get_axis_rz
B_STS_get_axis_rz.restype = a5err
B_STS_get_axis_rz.argtypes = [ctypes.c_double * 2, ctypes.POINTER(struct_c__SA_B_STS_data), real]
class struct_c__SA_B_TC_offload_data(Structure):
    pass

struct_c__SA_B_TC_offload_data._pack_ = 1 # source:False
struct_c__SA_B_TC_offload_data._fields_ = [
    ('axisr', ctypes.c_double),
    ('axisz', ctypes.c_double),
    ('psival', ctypes.c_double),
    ('rhoval', ctypes.c_double),
    ('B', ctypes.c_double * 3),
    ('dB', ctypes.c_double * 9),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
]

B_TC_offload_data = struct_c__SA_B_TC_offload_data
class struct_c__SA_B_TC_data(Structure):
    pass

struct_c__SA_B_TC_data._pack_ = 1 # source:False
struct_c__SA_B_TC_data._fields_ = [
    ('axisr', ctypes.c_double),
    ('axisz', ctypes.c_double),
    ('psival', ctypes.c_double),
    ('rhoval', ctypes.c_double),
    ('B', ctypes.POINTER(ctypes.c_double)),
    ('dB', ctypes.POINTER(ctypes.c_double)),
]

B_TC_data = struct_c__SA_B_TC_data
B_TC_init_offload = _libraries['libascot.so'].B_TC_init_offload
B_TC_init_offload.restype = ctypes.c_int32
B_TC_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_TC_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_TC_free_offload = _libraries['libascot.so'].B_TC_free_offload
B_TC_free_offload.restype = None
B_TC_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_TC_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_TC_init = _libraries['libascot.so'].B_TC_init
B_TC_init.restype = None
B_TC_init.argtypes = [ctypes.POINTER(struct_c__SA_B_TC_data), ctypes.POINTER(struct_c__SA_B_TC_offload_data), ctypes.POINTER(ctypes.c_double)]
B_TC_eval_B = _libraries['libascot.so'].B_TC_eval_B
B_TC_eval_B.restype = a5err
B_TC_eval_B.argtypes = [ctypes.c_double * 3, real, real, real, ctypes.POINTER(struct_c__SA_B_TC_data)]
B_TC_eval_psi = _libraries['libascot.so'].B_TC_eval_psi
B_TC_eval_psi.restype = a5err
B_TC_eval_psi.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, ctypes.POINTER(struct_c__SA_B_TC_data)]
B_TC_eval_psi_dpsi = _libraries['libascot.so'].B_TC_eval_psi_dpsi
B_TC_eval_psi_dpsi.restype = a5err
B_TC_eval_psi_dpsi.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_TC_data)]
B_TC_eval_rho_drho = _libraries['libascot.so'].B_TC_eval_rho_drho
B_TC_eval_rho_drho.restype = a5err
B_TC_eval_rho_drho.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_TC_data)]
B_TC_eval_B_dB = _libraries['libascot.so'].B_TC_eval_B_dB
B_TC_eval_B_dB.restype = a5err
B_TC_eval_B_dB.argtypes = [ctypes.c_double * 12, real, real, real, ctypes.POINTER(struct_c__SA_B_TC_data)]
B_TC_get_axis_rz = _libraries['libascot.so'].B_TC_get_axis_rz
B_TC_get_axis_rz.restype = a5err
B_TC_get_axis_rz.argtypes = [ctypes.c_double * 2, ctypes.POINTER(struct_c__SA_B_TC_data)]

# values for enumeration 'B_field_type'
B_field_type__enumvalues = {
    0: 'B_field_type_GS',
    1: 'B_field_type_2DS',
    2: 'B_field_type_3DS',
    3: 'B_field_type_STS',
    4: 'B_field_type_TC',
}
B_field_type_GS = 0
B_field_type_2DS = 1
B_field_type_3DS = 2
B_field_type_STS = 3
B_field_type_TC = 4
B_field_type = ctypes.c_uint32 # enum
class struct_c__SA_B_field_offload_data(Structure):
    pass

struct_c__SA_B_field_offload_data._pack_ = 1 # source:False
struct_c__SA_B_field_offload_data._fields_ = [
    ('type', B_field_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('BGS', B_GS_offload_data),
    ('B2DS', B_2DS_offload_data),
    ('B3DS', B_3DS_offload_data),
    ('BSTS', B_STS_offload_data),
    ('BTC', B_TC_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

B_field_offload_data = struct_c__SA_B_field_offload_data
class struct_c__SA_B_field_data(Structure):
    pass

struct_c__SA_B_field_data._pack_ = 1 # source:False
struct_c__SA_B_field_data._fields_ = [
    ('type', B_field_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('BGS', B_GS_data),
    ('B2DS', B_2DS_data),
    ('B3DS', B_3DS_data),
    ('BSTS', B_STS_data),
    ('BTC', B_TC_data),
]

B_field_data = struct_c__SA_B_field_data
B_field_init_offload = _libraries['libascot.so'].B_field_init_offload
B_field_init_offload.restype = ctypes.c_int32
B_field_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_field_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_field_free_offload = _libraries['libascot.so'].B_field_free_offload
B_field_free_offload.restype = None
B_field_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_B_field_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
B_field_init = _libraries['libascot.so'].B_field_init
B_field_init.restype = ctypes.c_int32
B_field_init.argtypes = [ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(struct_c__SA_B_field_offload_data), ctypes.POINTER(ctypes.c_double)]
B_field_eval_psi = _libraries['libascot.so'].B_field_eval_psi
B_field_eval_psi.restype = a5err
B_field_eval_psi.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, real, ctypes.POINTER(struct_c__SA_B_field_data)]
B_field_eval_psi_dpsi = _libraries['libascot.so'].B_field_eval_psi_dpsi
B_field_eval_psi_dpsi.restype = a5err
B_field_eval_psi_dpsi.argtypes = [ctypes.c_double * 4, real, real, real, real, ctypes.POINTER(struct_c__SA_B_field_data)]
B_field_eval_rho = _libraries['libascot.so'].B_field_eval_rho
B_field_eval_rho.restype = a5err
B_field_eval_rho.argtypes = [ctypes.c_double * 2, real, ctypes.POINTER(struct_c__SA_B_field_data)]
B_field_eval_rho_drho = _libraries['libascot.so'].B_field_eval_rho_drho
B_field_eval_rho_drho.restype = a5err
B_field_eval_rho_drho.argtypes = [ctypes.c_double * 4, real, real, real, ctypes.POINTER(struct_c__SA_B_field_data)]
B_field_eval_B = _libraries['libascot.so'].B_field_eval_B
B_field_eval_B.restype = a5err
B_field_eval_B.argtypes = [ctypes.c_double * 3, real, real, real, real, ctypes.POINTER(struct_c__SA_B_field_data)]
B_field_eval_B_dB = _libraries['libascot.so'].B_field_eval_B_dB
B_field_eval_B_dB.restype = a5err
B_field_eval_B_dB.argtypes = [ctypes.c_double * 15, real, real, real, real, ctypes.POINTER(struct_c__SA_B_field_data)]
B_field_get_axis_rz = _libraries['libascot.so'].B_field_get_axis_rz
B_field_get_axis_rz.restype = a5err
B_field_get_axis_rz.argtypes = [ctypes.c_double * 2, ctypes.POINTER(struct_c__SA_B_field_data), real]
class struct_c__SA_E_TC_offload_data(Structure):
    pass

struct_c__SA_E_TC_offload_data._pack_ = 1 # source:False
struct_c__SA_E_TC_offload_data._fields_ = [
    ('Exyz', ctypes.c_double * 3),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
]

E_TC_offload_data = struct_c__SA_E_TC_offload_data
class struct_c__SA_E_TC_data(Structure):
    pass

struct_c__SA_E_TC_data._pack_ = 1 # source:False
struct_c__SA_E_TC_data._fields_ = [
    ('Exyz', ctypes.POINTER(ctypes.c_double)),
]

E_TC_data = struct_c__SA_E_TC_data
E_TC_init_offload = _libraries['libascot.so'].E_TC_init_offload
E_TC_init_offload.restype = ctypes.c_int32
E_TC_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_E_TC_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
E_TC_free_offload = _libraries['libascot.so'].E_TC_free_offload
E_TC_free_offload.restype = None
E_TC_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_E_TC_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
E_TC_init = _libraries['libascot.so'].E_TC_init
E_TC_init.restype = None
E_TC_init.argtypes = [ctypes.POINTER(struct_c__SA_E_TC_data), ctypes.POINTER(struct_c__SA_E_TC_offload_data), ctypes.POINTER(ctypes.c_double)]
E_TC_eval_E = _libraries['libascot.so'].E_TC_eval_E
E_TC_eval_E.restype = a5err
E_TC_eval_E.argtypes = [ctypes.c_double * 3, real, real, real, ctypes.POINTER(struct_c__SA_E_TC_data), ctypes.POINTER(struct_c__SA_B_field_data)]
class struct_c__SA_E_1DS_offload_data(Structure):
    pass

struct_c__SA_E_1DS_offload_data._pack_ = 1 # source:False
struct_c__SA_E_1DS_offload_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('rho_min', ctypes.c_double),
    ('rho_max', ctypes.c_double),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

E_1DS_offload_data = struct_c__SA_E_1DS_offload_data
class struct_c__SA_E_1DS_data(Structure):
    pass

class struct_c__SA_interp1D_data(Structure):
    pass

struct_c__SA_interp1D_data._pack_ = 1 # source:False
struct_c__SA_interp1D_data._fields_ = [
    ('n_x', ctypes.c_int32),
    ('bc_x', ctypes.c_int32),
    ('x_min', ctypes.c_double),
    ('x_max', ctypes.c_double),
    ('x_grid', ctypes.c_double),
    ('c', ctypes.POINTER(ctypes.c_double)),
]

struct_c__SA_E_1DS_data._pack_ = 1 # source:False
struct_c__SA_E_1DS_data._fields_ = [
    ('dV', struct_c__SA_interp1D_data),
]

E_1DS_data = struct_c__SA_E_1DS_data
E_1DS_init_offload = _libraries['libascot.so'].E_1DS_init_offload
E_1DS_init_offload.restype = ctypes.c_int32
E_1DS_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_E_1DS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
E_1DS_free_offload = _libraries['libascot.so'].E_1DS_free_offload
E_1DS_free_offload.restype = None
E_1DS_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_E_1DS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
E_1DS_init = _libraries['libascot.so'].E_1DS_init
E_1DS_init.restype = None
E_1DS_init.argtypes = [ctypes.POINTER(struct_c__SA_E_1DS_data), ctypes.POINTER(struct_c__SA_E_1DS_offload_data), ctypes.POINTER(ctypes.c_double)]
E_1DS_eval_E = _libraries['libascot.so'].E_1DS_eval_E
E_1DS_eval_E.restype = a5err
E_1DS_eval_E.argtypes = [ctypes.c_double * 3, real, real, real, ctypes.POINTER(struct_c__SA_E_1DS_data), ctypes.POINTER(struct_c__SA_B_field_data)]

# values for enumeration 'E_field_type'
E_field_type__enumvalues = {
    0: 'E_field_type_TC',
    1: 'E_field_type_1DS',
}
E_field_type_TC = 0
E_field_type_1DS = 1
E_field_type = ctypes.c_uint32 # enum
class struct_c__SA_E_field_offload_data(Structure):
    pass

struct_c__SA_E_field_offload_data._pack_ = 1 # source:False
struct_c__SA_E_field_offload_data._fields_ = [
    ('type', E_field_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('ETC', E_TC_offload_data),
    ('E1DS', E_1DS_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

E_field_offload_data = struct_c__SA_E_field_offload_data
class struct_c__SA_E_field_data(Structure):
    pass

struct_c__SA_E_field_data._pack_ = 1 # source:False
struct_c__SA_E_field_data._fields_ = [
    ('type', E_field_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('ETC', E_TC_data),
    ('E1DS', E_1DS_data),
]

E_field_data = struct_c__SA_E_field_data
E_field_init_offload = _libraries['libascot.so'].E_field_init_offload
E_field_init_offload.restype = ctypes.c_int32
E_field_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_E_field_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
E_field_free_offload = _libraries['libascot.so'].E_field_free_offload
E_field_free_offload.restype = None
E_field_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_E_field_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
E_field_init = _libraries['libascot.so'].E_field_init
E_field_init.restype = ctypes.c_int32
E_field_init.argtypes = [ctypes.POINTER(struct_c__SA_E_field_data), ctypes.POINTER(struct_c__SA_E_field_offload_data), ctypes.POINTER(ctypes.c_double)]
E_field_eval_E = _libraries['libascot.so'].E_field_eval_E
E_field_eval_E.restype = a5err
E_field_eval_E.argtypes = [ctypes.c_double * 3, real, real, real, real, ctypes.POINTER(struct_c__SA_E_field_data), ctypes.POINTER(struct_c__SA_B_field_data)]
class struct_c__SA_particle_state(Structure):
    pass

struct_c__SA_particle_state._pack_ = 1 # source:False
struct_c__SA_particle_state._fields_ = [
    ('r', ctypes.c_double),
    ('phi', ctypes.c_double),
    ('z', ctypes.c_double),
    ('ppar', ctypes.c_double),
    ('mu', ctypes.c_double),
    ('zeta', ctypes.c_double),
    ('rprt', ctypes.c_double),
    ('phiprt', ctypes.c_double),
    ('zprt', ctypes.c_double),
    ('p_r', ctypes.c_double),
    ('p_phi', ctypes.c_double),
    ('p_z', ctypes.c_double),
    ('mass', ctypes.c_double),
    ('charge', ctypes.c_double),
    ('anum', ctypes.c_int32),
    ('znum', ctypes.c_int32),
    ('weight', ctypes.c_double),
    ('time', ctypes.c_double),
    ('mileage', ctypes.c_double),
    ('cputime', ctypes.c_double),
    ('rho', ctypes.c_double),
    ('theta', ctypes.c_double),
    ('id', ctypes.c_int64),
    ('endcond', ctypes.c_int64),
    ('walltile', ctypes.c_int64),
    ('B_r', ctypes.c_double),
    ('B_phi', ctypes.c_double),
    ('B_z', ctypes.c_double),
    ('B_r_dr', ctypes.c_double),
    ('B_phi_dr', ctypes.c_double),
    ('B_z_dr', ctypes.c_double),
    ('B_r_dphi', ctypes.c_double),
    ('B_phi_dphi', ctypes.c_double),
    ('B_z_dphi', ctypes.c_double),
    ('B_r_dz', ctypes.c_double),
    ('B_phi_dz', ctypes.c_double),
    ('B_z_dz', ctypes.c_double),
    ('err', ctypes.c_uint64),
]

particle_state = struct_c__SA_particle_state
class struct_c__SA_particle(Structure):
    pass

struct_c__SA_particle._pack_ = 1 # source:False
struct_c__SA_particle._fields_ = [
    ('r', ctypes.c_double),
    ('phi', ctypes.c_double),
    ('z', ctypes.c_double),
    ('p_r', ctypes.c_double),
    ('p_phi', ctypes.c_double),
    ('p_z', ctypes.c_double),
    ('mass', ctypes.c_double),
    ('charge', ctypes.c_double),
    ('anum', ctypes.c_int32),
    ('znum', ctypes.c_int32),
    ('weight', ctypes.c_double),
    ('time', ctypes.c_double),
    ('id', ctypes.c_int64),
]

particle = struct_c__SA_particle
class struct_c__SA_particle_gc(Structure):
    pass

struct_c__SA_particle_gc._pack_ = 1 # source:False
struct_c__SA_particle_gc._fields_ = [
    ('r', ctypes.c_double),
    ('phi', ctypes.c_double),
    ('z', ctypes.c_double),
    ('energy', ctypes.c_double),
    ('pitch', ctypes.c_double),
    ('zeta', ctypes.c_double),
    ('mass', ctypes.c_double),
    ('charge', ctypes.c_double),
    ('anum', ctypes.c_int32),
    ('znum', ctypes.c_int32),
    ('weight', ctypes.c_double),
    ('time', ctypes.c_double),
    ('id', ctypes.c_int64),
]

particle_gc = struct_c__SA_particle_gc
class struct_c__SA_particle_ml(Structure):
    pass

struct_c__SA_particle_ml._pack_ = 1 # source:False
struct_c__SA_particle_ml._fields_ = [
    ('r', ctypes.c_double),
    ('phi', ctypes.c_double),
    ('z', ctypes.c_double),
    ('pitch', ctypes.c_double),
    ('weight', ctypes.c_double),
    ('time', ctypes.c_double),
    ('id', ctypes.c_int64),
]

particle_ml = struct_c__SA_particle_ml
class struct_c__SA_particle_queue(Structure):
    pass

struct_c__SA_particle_queue._pack_ = 1 # source:False
struct_c__SA_particle_queue._fields_ = [
    ('n', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('p', ctypes.POINTER(ctypes.POINTER(struct_c__SA_particle_state))),
    ('next', ctypes.c_int32),
    ('finished', ctypes.c_int32),
]

particle_queue = struct_c__SA_particle_queue

# values for enumeration 'input_particle_type'
input_particle_type__enumvalues = {
    0: 'input_particle_type_p',
    1: 'input_particle_type_gc',
    2: 'input_particle_type_ml',
    3: 'input_particle_type_s',
}
input_particle_type_p = 0
input_particle_type_gc = 1
input_particle_type_ml = 2
input_particle_type_s = 3
input_particle_type = ctypes.c_uint32 # enum
class struct_c__SA_input_particle(Structure):
    pass

class union_c__SA_input_particle_0(Union):
    _pack_ = 1 # source:False
    _fields_ = [
    ('p', particle),
    ('p_gc', particle_gc),
    ('p_ml', particle_ml),
    ('p_s', particle_state),
     ]

struct_c__SA_input_particle._pack_ = 1 # source:False
struct_c__SA_input_particle._anonymous_ = ('_0',)
struct_c__SA_input_particle._fields_ = [
    ('type', input_particle_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('_0', union_c__SA_input_particle_0),
]

input_particle = struct_c__SA_input_particle
class struct_c__SA_particle_simd_fo(Structure):
    pass

struct_c__SA_particle_simd_fo._pack_ = 1 # source:False
struct_c__SA_particle_simd_fo._fields_ = [
    ('r', ctypes.c_double * 16),
    ('phi', ctypes.c_double * 16),
    ('z', ctypes.c_double * 16),
    ('p_r', ctypes.c_double * 16),
    ('p_phi', ctypes.c_double * 16),
    ('p_z', ctypes.c_double * 16),
    ('mass', ctypes.c_double * 16),
    ('charge', ctypes.c_double * 16),
    ('time', ctypes.c_double * 16),
    ('znum', ctypes.c_int32 * 16),
    ('anum', ctypes.c_int32 * 16),
    ('B_r', ctypes.c_double * 16),
    ('B_phi', ctypes.c_double * 16),
    ('B_z', ctypes.c_double * 16),
    ('B_r_dr', ctypes.c_double * 16),
    ('B_phi_dr', ctypes.c_double * 16),
    ('B_z_dr', ctypes.c_double * 16),
    ('B_r_dphi', ctypes.c_double * 16),
    ('B_phi_dphi', ctypes.c_double * 16),
    ('B_z_dphi', ctypes.c_double * 16),
    ('B_r_dz', ctypes.c_double * 16),
    ('B_phi_dz', ctypes.c_double * 16),
    ('B_z_dz', ctypes.c_double * 16),
    ('bounces', ctypes.c_int32 * 16),
    ('weight', ctypes.c_double * 16),
    ('cputime', ctypes.c_double * 16),
    ('rho', ctypes.c_double * 16),
    ('theta', ctypes.c_double * 16),
    ('id', ctypes.c_int64 * 16),
    ('endcond', ctypes.c_int64 * 16),
    ('walltile', ctypes.c_int64 * 16),
    ('mileage', ctypes.c_double * 16),
    ('running', ctypes.c_int64 * 16),
    ('err', ctypes.c_uint64 * 16),
    ('index', ctypes.c_int64 * 16),
]

particle_simd_fo = struct_c__SA_particle_simd_fo
class struct_c__SA_particle_simd_gc(Structure):
    pass

struct_c__SA_particle_simd_gc._pack_ = 1 # source:False
struct_c__SA_particle_simd_gc._fields_ = [
    ('r', ctypes.c_double * 16),
    ('phi', ctypes.c_double * 16),
    ('z', ctypes.c_double * 16),
    ('ppar', ctypes.c_double * 16),
    ('mu', ctypes.c_double * 16),
    ('zeta', ctypes.c_double * 16),
    ('mass', ctypes.c_double * 16),
    ('charge', ctypes.c_double * 16),
    ('time', ctypes.c_double * 16),
    ('B_r', ctypes.c_double * 16),
    ('B_phi', ctypes.c_double * 16),
    ('B_z', ctypes.c_double * 16),
    ('B_r_dr', ctypes.c_double * 16),
    ('B_phi_dr', ctypes.c_double * 16),
    ('B_z_dr', ctypes.c_double * 16),
    ('B_r_dphi', ctypes.c_double * 16),
    ('B_phi_dphi', ctypes.c_double * 16),
    ('B_z_dphi', ctypes.c_double * 16),
    ('B_r_dz', ctypes.c_double * 16),
    ('B_phi_dz', ctypes.c_double * 16),
    ('B_z_dz', ctypes.c_double * 16),
    ('bounces', ctypes.c_int32 * 16),
    ('weight', ctypes.c_double * 16),
    ('cputime', ctypes.c_double * 16),
    ('rho', ctypes.c_double * 16),
    ('theta', ctypes.c_double * 16),
    ('id', ctypes.c_int64 * 16),
    ('endcond', ctypes.c_int64 * 16),
    ('walltile', ctypes.c_int64 * 16),
    ('mileage', ctypes.c_double * 16),
    ('running', ctypes.c_int64 * 16),
    ('err', ctypes.c_uint64 * 16),
    ('index', ctypes.c_int64 * 16),
]

particle_simd_gc = struct_c__SA_particle_simd_gc
class struct_c__SA_particle_simd_ml(Structure):
    pass

struct_c__SA_particle_simd_ml._pack_ = 1 # source:False
struct_c__SA_particle_simd_ml._fields_ = [
    ('r', ctypes.c_double * 16),
    ('phi', ctypes.c_double * 16),
    ('z', ctypes.c_double * 16),
    ('pitch', ctypes.c_double * 16),
    ('time', ctypes.c_double * 16),
    ('B_r', ctypes.c_double * 16),
    ('B_phi', ctypes.c_double * 16),
    ('B_z', ctypes.c_double * 16),
    ('B_r_dr', ctypes.c_double * 16),
    ('B_phi_dr', ctypes.c_double * 16),
    ('B_z_dr', ctypes.c_double * 16),
    ('B_r_dphi', ctypes.c_double * 16),
    ('B_phi_dphi', ctypes.c_double * 16),
    ('B_z_dphi', ctypes.c_double * 16),
    ('B_r_dz', ctypes.c_double * 16),
    ('B_phi_dz', ctypes.c_double * 16),
    ('B_z_dz', ctypes.c_double * 16),
    ('weight', ctypes.c_double * 16),
    ('cputime', ctypes.c_double * 16),
    ('rho', ctypes.c_double * 16),
    ('theta', ctypes.c_double * 16),
    ('id', ctypes.c_int64 * 16),
    ('endcond', ctypes.c_int64 * 16),
    ('walltile', ctypes.c_int64 * 16),
    ('mileage', ctypes.c_double * 16),
    ('running', ctypes.c_int64 * 16),
    ('err', ctypes.c_uint64 * 16),
    ('index', ctypes.c_int64 * 16),
]

particle_simd_ml = struct_c__SA_particle_simd_ml
particle_to_fo_dummy = _libraries['libascot.so'].particle_to_fo_dummy
particle_to_fo_dummy.restype = None
particle_to_fo_dummy.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.c_int32]
particle_to_gc_dummy = _libraries['libascot.so'].particle_to_gc_dummy
particle_to_gc_dummy.restype = None
particle_to_gc_dummy.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.c_int32]
particle_to_ml_dummy = _libraries['libascot.so'].particle_to_ml_dummy
particle_to_ml_dummy.restype = None
particle_to_ml_dummy.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.c_int32]
particle_cycle_fo = _libraries['libascot.so'].particle_cycle_fo
particle_cycle_fo.restype = ctypes.c_int32
particle_cycle_fo.argtypes = [ctypes.POINTER(struct_c__SA_particle_queue), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(ctypes.c_int32)]
particle_cycle_gc = _libraries['libascot.so'].particle_cycle_gc
particle_cycle_gc.restype = ctypes.c_int32
particle_cycle_gc.argtypes = [ctypes.POINTER(struct_c__SA_particle_queue), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(ctypes.c_int32)]
particle_cycle_ml = _libraries['libascot.so'].particle_cycle_ml
particle_cycle_ml.restype = ctypes.c_int32
particle_cycle_ml.argtypes = [ctypes.POINTER(struct_c__SA_particle_queue), ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(ctypes.c_int32)]
particle_input_to_state = _libraries['libascot.so'].particle_input_to_state
particle_input_to_state.restype = None
particle_input_to_state.argtypes = [ctypes.POINTER(struct_c__SA_input_particle), ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_B_field_data)]
particle_input_p_to_state = _libraries['libascot.so'].particle_input_p_to_state
particle_input_p_to_state.restype = a5err
particle_input_p_to_state.argtypes = [ctypes.POINTER(struct_c__SA_particle), ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_B_field_data)]
particle_input_gc_to_state = _libraries['libascot.so'].particle_input_gc_to_state
particle_input_gc_to_state.restype = a5err
particle_input_gc_to_state.argtypes = [ctypes.POINTER(struct_c__SA_particle_gc), ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_B_field_data)]
particle_input_ml_to_state = _libraries['libascot.so'].particle_input_ml_to_state
particle_input_ml_to_state.restype = a5err
particle_input_ml_to_state.argtypes = [ctypes.POINTER(struct_c__SA_particle_ml), ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_B_field_data)]
particle_state_to_fo = _libraries['libascot.so'].particle_state_to_fo
particle_state_to_fo.restype = a5err
particle_state_to_fo.argtypes = [ctypes.POINTER(struct_c__SA_particle_state), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.c_int32, ctypes.POINTER(struct_c__SA_B_field_data)]
particle_fo_to_state = _libraries['libascot.so'].particle_fo_to_state
particle_fo_to_state.restype = None
particle_fo_to_state.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_B_field_data)]
particle_state_to_gc = _libraries['libascot.so'].particle_state_to_gc
particle_state_to_gc.restype = a5err
particle_state_to_gc.argtypes = [ctypes.POINTER(struct_c__SA_particle_state), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.c_int32, ctypes.POINTER(struct_c__SA_B_field_data)]
particle_gc_to_state = _libraries['libascot.so'].particle_gc_to_state
particle_gc_to_state.restype = None
particle_gc_to_state.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_B_field_data)]
particle_state_to_ml = _libraries['libascot.so'].particle_state_to_ml
particle_state_to_ml.restype = a5err
particle_state_to_ml.argtypes = [ctypes.POINTER(struct_c__SA_particle_state), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.c_int32, ctypes.POINTER(struct_c__SA_B_field_data)]
particle_ml_to_state = _libraries['libascot.so'].particle_ml_to_state
particle_ml_to_state.restype = None
particle_ml_to_state.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_B_field_data)]
particle_fo_to_gc = _libraries['libascot.so'].particle_fo_to_gc
particle_fo_to_gc.restype = ctypes.c_int32
particle_fo_to_gc.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_B_field_data)]
particle_copy_fo = _libraries['libascot.so'].particle_copy_fo
particle_copy_fo.restype = None
particle_copy_fo.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.c_int32]
particle_copy_gc = _libraries['libascot.so'].particle_copy_gc
particle_copy_gc.restype = None
particle_copy_gc.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.c_int32]
particle_copy_ml = _libraries['libascot.so'].particle_copy_ml
particle_copy_ml.restype = None
particle_copy_ml.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.c_int32]
class struct_c__SA_dist_5D_offload_data(Structure):
    pass

struct_c__SA_dist_5D_offload_data._pack_ = 1 # source:False
struct_c__SA_dist_5D_offload_data._fields_ = [
    ('n_r', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_r', ctypes.c_double),
    ('max_r', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_z', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_z', ctypes.c_double),
    ('max_z', ctypes.c_double),
    ('n_ppara', ctypes.c_int32),
    ('PADDING_3', ctypes.c_ubyte * 4),
    ('min_ppara', ctypes.c_double),
    ('max_ppara', ctypes.c_double),
    ('n_pperp', ctypes.c_int32),
    ('PADDING_4', ctypes.c_ubyte * 4),
    ('min_pperp', ctypes.c_double),
    ('max_pperp', ctypes.c_double),
    ('n_time', ctypes.c_int32),
    ('PADDING_5', ctypes.c_ubyte * 4),
    ('min_time', ctypes.c_double),
    ('max_time', ctypes.c_double),
    ('n_q', ctypes.c_int32),
    ('PADDING_6', ctypes.c_ubyte * 4),
    ('min_q', ctypes.c_double),
    ('max_q', ctypes.c_double),
]

dist_5D_offload_data = struct_c__SA_dist_5D_offload_data
class struct_c__SA_dist_5D_data(Structure):
    pass

struct_c__SA_dist_5D_data._pack_ = 1 # source:False
struct_c__SA_dist_5D_data._fields_ = [
    ('n_r', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_r', ctypes.c_double),
    ('max_r', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_z', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_z', ctypes.c_double),
    ('max_z', ctypes.c_double),
    ('n_ppara', ctypes.c_int32),
    ('PADDING_3', ctypes.c_ubyte * 4),
    ('min_ppara', ctypes.c_double),
    ('max_ppara', ctypes.c_double),
    ('n_pperp', ctypes.c_int32),
    ('PADDING_4', ctypes.c_ubyte * 4),
    ('min_pperp', ctypes.c_double),
    ('max_pperp', ctypes.c_double),
    ('n_time', ctypes.c_int32),
    ('PADDING_5', ctypes.c_ubyte * 4),
    ('min_time', ctypes.c_double),
    ('max_time', ctypes.c_double),
    ('n_q', ctypes.c_int32),
    ('PADDING_6', ctypes.c_ubyte * 4),
    ('min_q', ctypes.c_double),
    ('max_q', ctypes.c_double),
    ('step_1', ctypes.c_uint64),
    ('step_2', ctypes.c_uint64),
    ('step_3', ctypes.c_uint64),
    ('step_4', ctypes.c_uint64),
    ('step_5', ctypes.c_uint64),
    ('step_6', ctypes.c_uint64),
    ('histogram', ctypes.POINTER(ctypes.c_double)),
]

dist_5D_data = struct_c__SA_dist_5D_data
size_t = ctypes.c_uint64
dist_5D_index = _libraries['libascot.so'].dist_5D_index
dist_5D_index.restype = size_t
dist_5D_index.argtypes = [ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, size_t, size_t, size_t, size_t, size_t, size_t]
dist_5D_init = _libraries['libascot.so'].dist_5D_init
dist_5D_init.restype = None
dist_5D_init.argtypes = [ctypes.POINTER(struct_c__SA_dist_5D_data), ctypes.POINTER(struct_c__SA_dist_5D_offload_data), ctypes.POINTER(ctypes.c_double)]
dist_5D_update_fo = _libraries['libascot.so'].dist_5D_update_fo
dist_5D_update_fo.restype = None
dist_5D_update_fo.argtypes = [ctypes.POINTER(struct_c__SA_dist_5D_data), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo)]
dist_5D_update_gc = _libraries['libascot.so'].dist_5D_update_gc
dist_5D_update_gc.restype = None
dist_5D_update_gc.argtypes = [ctypes.POINTER(struct_c__SA_dist_5D_data), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc)]
class struct_c__SA_dist_6D_offload_data(Structure):
    pass

struct_c__SA_dist_6D_offload_data._pack_ = 1 # source:False
struct_c__SA_dist_6D_offload_data._fields_ = [
    ('n_r', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_r', ctypes.c_double),
    ('max_r', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_z', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_z', ctypes.c_double),
    ('max_z', ctypes.c_double),
    ('n_pr', ctypes.c_int32),
    ('PADDING_3', ctypes.c_ubyte * 4),
    ('min_pr', ctypes.c_double),
    ('max_pr', ctypes.c_double),
    ('n_pphi', ctypes.c_int32),
    ('PADDING_4', ctypes.c_ubyte * 4),
    ('min_pphi', ctypes.c_double),
    ('max_pphi', ctypes.c_double),
    ('n_pz', ctypes.c_int32),
    ('PADDING_5', ctypes.c_ubyte * 4),
    ('min_pz', ctypes.c_double),
    ('max_pz', ctypes.c_double),
    ('n_time', ctypes.c_int32),
    ('PADDING_6', ctypes.c_ubyte * 4),
    ('min_time', ctypes.c_double),
    ('max_time', ctypes.c_double),
    ('n_q', ctypes.c_int32),
    ('PADDING_7', ctypes.c_ubyte * 4),
    ('min_q', ctypes.c_double),
    ('max_q', ctypes.c_double),
]

dist_6D_offload_data = struct_c__SA_dist_6D_offload_data
class struct_c__SA_dist_6D_data(Structure):
    pass

struct_c__SA_dist_6D_data._pack_ = 1 # source:False
struct_c__SA_dist_6D_data._fields_ = [
    ('n_r', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_r', ctypes.c_double),
    ('max_r', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_z', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_z', ctypes.c_double),
    ('max_z', ctypes.c_double),
    ('n_pr', ctypes.c_int32),
    ('PADDING_3', ctypes.c_ubyte * 4),
    ('min_pr', ctypes.c_double),
    ('max_pr', ctypes.c_double),
    ('n_pphi', ctypes.c_int32),
    ('PADDING_4', ctypes.c_ubyte * 4),
    ('min_pphi', ctypes.c_double),
    ('max_pphi', ctypes.c_double),
    ('n_pz', ctypes.c_int32),
    ('PADDING_5', ctypes.c_ubyte * 4),
    ('min_pz', ctypes.c_double),
    ('max_pz', ctypes.c_double),
    ('n_time', ctypes.c_int32),
    ('PADDING_6', ctypes.c_ubyte * 4),
    ('min_time', ctypes.c_double),
    ('max_time', ctypes.c_double),
    ('n_q', ctypes.c_int32),
    ('PADDING_7', ctypes.c_ubyte * 4),
    ('min_q', ctypes.c_double),
    ('max_q', ctypes.c_double),
    ('step_1', ctypes.c_uint64),
    ('step_2', ctypes.c_uint64),
    ('step_3', ctypes.c_uint64),
    ('step_4', ctypes.c_uint64),
    ('step_5', ctypes.c_uint64),
    ('step_6', ctypes.c_uint64),
    ('step_7', ctypes.c_uint64),
    ('histogram', ctypes.POINTER(ctypes.c_double)),
]

dist_6D_data = struct_c__SA_dist_6D_data
dist_6D_init = _libraries['libascot.so'].dist_6D_init
dist_6D_init.restype = None
dist_6D_init.argtypes = [ctypes.POINTER(struct_c__SA_dist_6D_data), ctypes.POINTER(struct_c__SA_dist_6D_offload_data), ctypes.POINTER(ctypes.c_double)]
dist_6D_update_fo = _libraries['libascot.so'].dist_6D_update_fo
dist_6D_update_fo.restype = None
dist_6D_update_fo.argtypes = [ctypes.POINTER(struct_c__SA_dist_6D_data), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo)]
dist_6D_update_gc = _libraries['libascot.so'].dist_6D_update_gc
dist_6D_update_gc.restype = None
dist_6D_update_gc.argtypes = [ctypes.POINTER(struct_c__SA_dist_6D_data), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc)]
class struct_c__SA_dist_rho5D_offload_data(Structure):
    pass

struct_c__SA_dist_rho5D_offload_data._pack_ = 1 # source:False
struct_c__SA_dist_rho5D_offload_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_rho', ctypes.c_double),
    ('max_rho', ctypes.c_double),
    ('n_theta', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_theta', ctypes.c_double),
    ('max_theta', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_ppara', ctypes.c_int32),
    ('PADDING_3', ctypes.c_ubyte * 4),
    ('min_ppara', ctypes.c_double),
    ('max_ppara', ctypes.c_double),
    ('n_pperp', ctypes.c_int32),
    ('PADDING_4', ctypes.c_ubyte * 4),
    ('min_pperp', ctypes.c_double),
    ('max_pperp', ctypes.c_double),
    ('n_time', ctypes.c_int32),
    ('PADDING_5', ctypes.c_ubyte * 4),
    ('min_time', ctypes.c_double),
    ('max_time', ctypes.c_double),
    ('n_q', ctypes.c_int32),
    ('PADDING_6', ctypes.c_ubyte * 4),
    ('min_q', ctypes.c_double),
    ('max_q', ctypes.c_double),
]

dist_rho5D_offload_data = struct_c__SA_dist_rho5D_offload_data
class struct_c__SA_dist_rho5D_data(Structure):
    pass

struct_c__SA_dist_rho5D_data._pack_ = 1 # source:False
struct_c__SA_dist_rho5D_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_rho', ctypes.c_double),
    ('max_rho', ctypes.c_double),
    ('n_theta', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_theta', ctypes.c_double),
    ('max_theta', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_ppara', ctypes.c_int32),
    ('PADDING_3', ctypes.c_ubyte * 4),
    ('min_ppara', ctypes.c_double),
    ('max_ppara', ctypes.c_double),
    ('n_pperp', ctypes.c_int32),
    ('PADDING_4', ctypes.c_ubyte * 4),
    ('min_pperp', ctypes.c_double),
    ('max_pperp', ctypes.c_double),
    ('n_time', ctypes.c_int32),
    ('PADDING_5', ctypes.c_ubyte * 4),
    ('min_time', ctypes.c_double),
    ('max_time', ctypes.c_double),
    ('n_q', ctypes.c_int32),
    ('PADDING_6', ctypes.c_ubyte * 4),
    ('min_q', ctypes.c_double),
    ('max_q', ctypes.c_double),
    ('step_1', ctypes.c_uint64),
    ('step_2', ctypes.c_uint64),
    ('step_3', ctypes.c_uint64),
    ('step_4', ctypes.c_uint64),
    ('step_5', ctypes.c_uint64),
    ('step_6', ctypes.c_uint64),
    ('histogram', ctypes.POINTER(ctypes.c_double)),
]

dist_rho5D_data = struct_c__SA_dist_rho5D_data
dist_rho5D_init = _libraries['libascot.so'].dist_rho5D_init
dist_rho5D_init.restype = None
dist_rho5D_init.argtypes = [ctypes.POINTER(struct_c__SA_dist_rho5D_data), ctypes.POINTER(struct_c__SA_dist_rho5D_offload_data), ctypes.POINTER(ctypes.c_double)]
dist_rho5D_update_fo = _libraries['libascot.so'].dist_rho5D_update_fo
dist_rho5D_update_fo.restype = None
dist_rho5D_update_fo.argtypes = [ctypes.POINTER(struct_c__SA_dist_rho5D_data), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo)]
dist_rho5D_update_gc = _libraries['libascot.so'].dist_rho5D_update_gc
dist_rho5D_update_gc.restype = None
dist_rho5D_update_gc.argtypes = [ctypes.POINTER(struct_c__SA_dist_rho5D_data), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc)]
class struct_c__SA_dist_rho6D_offload_data(Structure):
    pass

struct_c__SA_dist_rho6D_offload_data._pack_ = 1 # source:False
struct_c__SA_dist_rho6D_offload_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_rho', ctypes.c_double),
    ('max_rho', ctypes.c_double),
    ('n_theta', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_theta', ctypes.c_double),
    ('max_theta', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_pr', ctypes.c_int32),
    ('PADDING_3', ctypes.c_ubyte * 4),
    ('min_pr', ctypes.c_double),
    ('max_pr', ctypes.c_double),
    ('n_pphi', ctypes.c_int32),
    ('PADDING_4', ctypes.c_ubyte * 4),
    ('min_pphi', ctypes.c_double),
    ('max_pphi', ctypes.c_double),
    ('n_pz', ctypes.c_int32),
    ('PADDING_5', ctypes.c_ubyte * 4),
    ('min_pz', ctypes.c_double),
    ('max_pz', ctypes.c_double),
    ('n_time', ctypes.c_int32),
    ('PADDING_6', ctypes.c_ubyte * 4),
    ('min_time', ctypes.c_double),
    ('max_time', ctypes.c_double),
    ('n_q', ctypes.c_int32),
    ('PADDING_7', ctypes.c_ubyte * 4),
    ('min_q', ctypes.c_double),
    ('max_q', ctypes.c_double),
]

dist_rho6D_offload_data = struct_c__SA_dist_rho6D_offload_data
class struct_c__SA_dist_rho6D_data(Structure):
    pass

struct_c__SA_dist_rho6D_data._pack_ = 1 # source:False
struct_c__SA_dist_rho6D_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_rho', ctypes.c_double),
    ('max_rho', ctypes.c_double),
    ('n_theta', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_theta', ctypes.c_double),
    ('max_theta', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_pr', ctypes.c_int32),
    ('PADDING_3', ctypes.c_ubyte * 4),
    ('min_pr', ctypes.c_double),
    ('max_pr', ctypes.c_double),
    ('n_pphi', ctypes.c_int32),
    ('PADDING_4', ctypes.c_ubyte * 4),
    ('min_pphi', ctypes.c_double),
    ('max_pphi', ctypes.c_double),
    ('n_pz', ctypes.c_int32),
    ('PADDING_5', ctypes.c_ubyte * 4),
    ('min_pz', ctypes.c_double),
    ('max_pz', ctypes.c_double),
    ('n_time', ctypes.c_int32),
    ('PADDING_6', ctypes.c_ubyte * 4),
    ('min_time', ctypes.c_double),
    ('max_time', ctypes.c_double),
    ('n_q', ctypes.c_int32),
    ('PADDING_7', ctypes.c_ubyte * 4),
    ('min_q', ctypes.c_double),
    ('max_q', ctypes.c_double),
    ('step_1', ctypes.c_uint64),
    ('step_2', ctypes.c_uint64),
    ('step_3', ctypes.c_uint64),
    ('step_4', ctypes.c_uint64),
    ('step_5', ctypes.c_uint64),
    ('step_6', ctypes.c_uint64),
    ('step_7', ctypes.c_uint64),
    ('histogram', ctypes.POINTER(ctypes.c_double)),
]

dist_rho6D_data = struct_c__SA_dist_rho6D_data
dist_rho6D_init = _libraries['libascot.so'].dist_rho6D_init
dist_rho6D_init.restype = None
dist_rho6D_init.argtypes = [ctypes.POINTER(struct_c__SA_dist_rho6D_data), ctypes.POINTER(struct_c__SA_dist_rho6D_offload_data), ctypes.POINTER(ctypes.c_double)]
dist_rho6D_update_fo = _libraries['libascot.so'].dist_rho6D_update_fo
dist_rho6D_update_fo.restype = None
dist_rho6D_update_fo.argtypes = [ctypes.POINTER(struct_c__SA_dist_rho6D_data), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo)]
dist_rho6D_update_gc = _libraries['libascot.so'].dist_rho6D_update_gc
dist_rho6D_update_gc.restype = None
dist_rho6D_update_gc.argtypes = [ctypes.POINTER(struct_c__SA_dist_rho6D_data), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc)]
class struct_c__SA_dist_COM_offload_data(Structure):
    pass

struct_c__SA_dist_COM_offload_data._pack_ = 1 # source:False
struct_c__SA_dist_COM_offload_data._fields_ = [
    ('n_mu', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_mu', ctypes.c_double),
    ('max_mu', ctypes.c_double),
    ('n_Ekin', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_Ekin', ctypes.c_double),
    ('max_Ekin', ctypes.c_double),
    ('n_Ptor', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_Ptor', ctypes.c_double),
    ('max_Ptor', ctypes.c_double),
]

dist_COM_offload_data = struct_c__SA_dist_COM_offload_data
class struct_c__SA_dist_COM_data(Structure):
    pass

struct_c__SA_dist_COM_data._pack_ = 1 # source:False
struct_c__SA_dist_COM_data._fields_ = [
    ('n_mu', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_mu', ctypes.c_double),
    ('max_mu', ctypes.c_double),
    ('n_Ekin', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_Ekin', ctypes.c_double),
    ('max_Ekin', ctypes.c_double),
    ('n_Ptor', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_Ptor', ctypes.c_double),
    ('max_Ptor', ctypes.c_double),
    ('step_1', ctypes.c_uint64),
    ('step_2', ctypes.c_uint64),
    ('histogram', ctypes.POINTER(ctypes.c_double)),
]

dist_COM_data = struct_c__SA_dist_COM_data
dist_COM_init = _libraries['libascot.so'].dist_COM_init
dist_COM_init.restype = None
dist_COM_init.argtypes = [ctypes.POINTER(struct_c__SA_dist_COM_data), ctypes.POINTER(struct_c__SA_dist_COM_offload_data), ctypes.POINTER(ctypes.c_double)]
dist_COM_update_fo = _libraries['libascot.so'].dist_COM_update_fo
dist_COM_update_fo.restype = None
dist_COM_update_fo.argtypes = [ctypes.POINTER(struct_c__SA_dist_COM_data), ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo)]
dist_COM_update_gc = _libraries['libascot.so'].dist_COM_update_gc
dist_COM_update_gc.restype = None
dist_COM_update_gc.argtypes = [ctypes.POINTER(struct_c__SA_dist_COM_data), ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc)]
diag_orb_check_plane_crossing = _libraries['libascot.so'].diag_orb_check_plane_crossing
diag_orb_check_plane_crossing.restype = real
diag_orb_check_plane_crossing.argtypes = [real, real, real]
diag_orb_check_radial_crossing = _libraries['libascot.so'].diag_orb_check_radial_crossing
diag_orb_check_radial_crossing.restype = real
diag_orb_check_radial_crossing.argtypes = [real, real, real]
class struct_c__SA_diag_orb_offload_data(Structure):
    pass

struct_c__SA_diag_orb_offload_data._pack_ = 1 # source:False
struct_c__SA_diag_orb_offload_data._fields_ = [
    ('record_mode', ctypes.c_int32),
    ('mode', ctypes.c_int32),
    ('Npnt', ctypes.c_int32),
    ('Nmrk', ctypes.c_int32),
    ('Nfld', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('writeInterval', ctypes.c_double),
    ('ntoroidalplots', ctypes.c_int32),
    ('npoloidalplots', ctypes.c_int32),
    ('nradialplots', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('toroidalangles', ctypes.c_double * 30),
    ('poloidalangles', ctypes.c_double * 30),
    ('radialdistances', ctypes.c_double * 30),
]

diag_orb_offload_data = struct_c__SA_diag_orb_offload_data
class struct_c__SA_diag_orb_data(Structure):
    pass

struct_c__SA_diag_orb_data._pack_ = 1 # source:False
struct_c__SA_diag_orb_data._fields_ = [
    ('id', ctypes.POINTER(ctypes.c_double)),
    ('mileage', ctypes.POINTER(ctypes.c_double)),
    ('r', ctypes.POINTER(ctypes.c_double)),
    ('phi', ctypes.POINTER(ctypes.c_double)),
    ('z', ctypes.POINTER(ctypes.c_double)),
    ('p_r', ctypes.POINTER(ctypes.c_double)),
    ('p_phi', ctypes.POINTER(ctypes.c_double)),
    ('p_z', ctypes.POINTER(ctypes.c_double)),
    ('ppar', ctypes.POINTER(ctypes.c_double)),
    ('mu', ctypes.POINTER(ctypes.c_double)),
    ('zeta', ctypes.POINTER(ctypes.c_double)),
    ('weight', ctypes.POINTER(ctypes.c_double)),
    ('charge', ctypes.POINTER(ctypes.c_double)),
    ('rho', ctypes.POINTER(ctypes.c_double)),
    ('theta', ctypes.POINTER(ctypes.c_double)),
    ('B_r', ctypes.POINTER(ctypes.c_double)),
    ('B_phi', ctypes.POINTER(ctypes.c_double)),
    ('B_z', ctypes.POINTER(ctypes.c_double)),
    ('simmode', ctypes.POINTER(ctypes.c_double)),
    ('pncrid', ctypes.POINTER(ctypes.c_double)),
    ('pncrdi', ctypes.POINTER(ctypes.c_double)),
    ('mrk_pnt', ctypes.POINTER(ctypes.c_int64)),
    ('mrk_recorded', ctypes.POINTER(ctypes.c_double)),
    ('mode', ctypes.c_int32),
    ('Npnt', ctypes.c_int32),
    ('Nmrk', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('writeInterval', ctypes.c_double),
    ('ntoroidalplots', ctypes.c_int32),
    ('npoloidalplots', ctypes.c_int32),
    ('nradialplots', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('toroidalangles', ctypes.c_double * 30),
    ('poloidalangles', ctypes.c_double * 30),
    ('radialdistances', ctypes.c_double * 30),
]

diag_orb_data = struct_c__SA_diag_orb_data
diag_orb_init = _libraries['libascot.so'].diag_orb_init
diag_orb_init.restype = None
diag_orb_init.argtypes = [ctypes.POINTER(struct_c__SA_diag_orb_data), ctypes.POINTER(struct_c__SA_diag_orb_offload_data), ctypes.POINTER(ctypes.c_double)]
diag_orb_free = _libraries['libascot.so'].diag_orb_free
diag_orb_free.restype = None
diag_orb_free.argtypes = [ctypes.POINTER(struct_c__SA_diag_orb_data)]
diag_orb_update_fo = _libraries['libascot.so'].diag_orb_update_fo
diag_orb_update_fo.restype = None
diag_orb_update_fo.argtypes = [ctypes.POINTER(struct_c__SA_diag_orb_data), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo)]
diag_orb_update_gc = _libraries['libascot.so'].diag_orb_update_gc
diag_orb_update_gc.restype = None
diag_orb_update_gc.argtypes = [ctypes.POINTER(struct_c__SA_diag_orb_data), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc)]
diag_orb_update_ml = _libraries['libascot.so'].diag_orb_update_ml
diag_orb_update_ml.restype = None
diag_orb_update_ml.argtypes = [ctypes.POINTER(struct_c__SA_diag_orb_data), ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.POINTER(struct_c__SA_particle_simd_ml)]
class struct_diag_transcoef_link(Structure):
    pass

struct_diag_transcoef_link._pack_ = 1 # source:False
struct_diag_transcoef_link._fields_ = [
    ('rho', ctypes.c_double),
    ('time', ctypes.c_double),
    ('pitchsign', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('prevlink', ctypes.POINTER(struct_diag_transcoef_link)),
]

diag_transcoef_link = struct_diag_transcoef_link
class struct_c__SA_diag_transcoef_offload_data(Structure):
    pass

struct_c__SA_diag_transcoef_offload_data._pack_ = 1 # source:False
struct_c__SA_diag_transcoef_offload_data._fields_ = [
    ('Nmrk', ctypes.c_int64),
    ('Navg', ctypes.c_int32),
    ('recordrho', ctypes.c_int32),
    ('interval', ctypes.c_double),
]

diag_transcoef_offload_data = struct_c__SA_diag_transcoef_offload_data
class struct_c__SA_diag_transcoef_data(Structure):
    pass

struct_c__SA_diag_transcoef_data._pack_ = 1 # source:False
struct_c__SA_diag_transcoef_data._fields_ = [
    ('Navg', ctypes.c_int32),
    ('recordrho', ctypes.c_int32),
    ('interval', ctypes.c_double),
    ('datapoints', ctypes.POINTER(ctypes.POINTER(struct_diag_transcoef_link))),
    ('id', ctypes.POINTER(ctypes.c_double)),
    ('Kcoef', ctypes.POINTER(ctypes.c_double)),
    ('Dcoef', ctypes.POINTER(ctypes.c_double)),
]

diag_transcoef_data = struct_c__SA_diag_transcoef_data
diag_transcoef_init = _libraries['libascot.so'].diag_transcoef_init
diag_transcoef_init.restype = None
diag_transcoef_init.argtypes = [ctypes.POINTER(struct_c__SA_diag_transcoef_data), ctypes.POINTER(struct_c__SA_diag_transcoef_offload_data), ctypes.POINTER(ctypes.c_double)]
diag_transcoef_free = _libraries['libascot.so'].diag_transcoef_free
diag_transcoef_free.restype = None
diag_transcoef_free.argtypes = [ctypes.POINTER(struct_c__SA_diag_transcoef_data)]
diag_transcoef_update_fo = _libraries['libascot.so'].diag_transcoef_update_fo
diag_transcoef_update_fo.restype = None
diag_transcoef_update_fo.argtypes = [ctypes.POINTER(struct_c__SA_diag_transcoef_data), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo)]
diag_transcoef_update_gc = _libraries['libascot.so'].diag_transcoef_update_gc
diag_transcoef_update_gc.restype = None
diag_transcoef_update_gc.argtypes = [ctypes.POINTER(struct_c__SA_diag_transcoef_data), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc)]
diag_transcoef_update_ml = _libraries['libascot.so'].diag_transcoef_update_ml
diag_transcoef_update_ml.restype = None
diag_transcoef_update_ml.argtypes = [ctypes.POINTER(struct_c__SA_diag_transcoef_data), ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.POINTER(struct_c__SA_particle_simd_ml)]
class struct_c__SA_diag_offload_data(Structure):
    pass

struct_c__SA_diag_offload_data._pack_ = 1 # source:False
struct_c__SA_diag_offload_data._fields_ = [
    ('diagorb_collect', ctypes.c_int32),
    ('dist5D_collect', ctypes.c_int32),
    ('dist6D_collect', ctypes.c_int32),
    ('distrho5D_collect', ctypes.c_int32),
    ('distrho6D_collect', ctypes.c_int32),
    ('distCOM_collect', ctypes.c_int32),
    ('diagtrcof_collect', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('diagorb', diag_orb_offload_data),
    ('dist5D', dist_5D_offload_data),
    ('dist6D', dist_6D_offload_data),
    ('distrho5D', dist_rho5D_offload_data),
    ('distrho6D', dist_rho6D_offload_data),
    ('distCOM', dist_COM_offload_data),
    ('diagtrcof', diag_transcoef_offload_data),
    ('offload_dist5D_index', ctypes.c_uint64),
    ('offload_dist6D_index', ctypes.c_uint64),
    ('offload_distrho5D_index', ctypes.c_uint64),
    ('offload_distrho6D_index', ctypes.c_uint64),
    ('offload_distCOM_index', ctypes.c_uint64),
    ('offload_diagorb_index', ctypes.c_uint64),
    ('offload_diagtrcof_index', ctypes.c_uint64),
    ('offload_dist_length', ctypes.c_uint64),
    ('offload_array_length', ctypes.c_uint64),
]

diag_offload_data = struct_c__SA_diag_offload_data
class struct_c__SA_diag_data(Structure):
    pass

struct_c__SA_diag_data._pack_ = 1 # source:False
struct_c__SA_diag_data._fields_ = [
    ('diagorb_collect', ctypes.c_int32),
    ('dist5D_collect', ctypes.c_int32),
    ('dist6D_collect', ctypes.c_int32),
    ('distrho5D_collect', ctypes.c_int32),
    ('distrho6D_collect', ctypes.c_int32),
    ('distCOM_collect', ctypes.c_int32),
    ('diagtrcof_collect', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('diagorb', diag_orb_data),
    ('dist5D', dist_5D_data),
    ('dist6D', dist_6D_data),
    ('distrho5D', dist_rho5D_data),
    ('distrho6D', dist_rho6D_data),
    ('distCOM', dist_COM_data),
    ('diagtrcof', diag_transcoef_data),
]

diag_data = struct_c__SA_diag_data
diag_init_offload = _libraries['libascot.so'].diag_init_offload
diag_init_offload.restype = ctypes.c_int32
diag_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_diag_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.c_int32]
diag_free_offload = _libraries['libascot.so'].diag_free_offload
diag_free_offload.restype = None
diag_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_diag_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
diag_sum = _libraries['libascot.so'].diag_sum
diag_sum.restype = None
diag_sum.argtypes = [ctypes.POINTER(struct_c__SA_diag_offload_data), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
diag_init = _libraries['libascot.so'].diag_init
diag_init.restype = None
diag_init.argtypes = [ctypes.POINTER(struct_c__SA_diag_data), ctypes.POINTER(struct_c__SA_diag_offload_data), ctypes.POINTER(ctypes.c_double)]
diag_free = _libraries['libascot.so'].diag_free
diag_free.restype = None
diag_free.argtypes = [ctypes.POINTER(struct_c__SA_diag_data)]
diag_update_fo = _libraries['libascot.so'].diag_update_fo
diag_update_fo.restype = None
diag_update_fo.argtypes = [ctypes.POINTER(struct_c__SA_diag_data), ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo)]
diag_update_gc = _libraries['libascot.so'].diag_update_gc
diag_update_gc.restype = None
diag_update_gc.argtypes = [ctypes.POINTER(struct_c__SA_diag_data), ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc)]
diag_update_ml = _libraries['libascot.so'].diag_update_ml
diag_update_ml.restype = None
diag_update_ml.argtypes = [ctypes.POINTER(struct_c__SA_diag_data), ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.POINTER(struct_c__SA_particle_simd_ml)]
class struct_c__SA_plasma_1D_offload_data(Structure):
    pass

struct_c__SA_plasma_1D_offload_data._pack_ = 1 # source:False
struct_c__SA_plasma_1D_offload_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('n_species', ctypes.c_int32),
    ('mass', ctypes.c_double * 8),
    ('charge', ctypes.c_double * 8),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
]

plasma_1D_offload_data = struct_c__SA_plasma_1D_offload_data
class struct_c__SA_plasma_1D_data(Structure):
    pass

struct_c__SA_plasma_1D_data._pack_ = 1 # source:False
struct_c__SA_plasma_1D_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('n_species', ctypes.c_int32),
    ('mass', ctypes.c_double * 8),
    ('charge', ctypes.c_double * 8),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('rho', ctypes.POINTER(ctypes.c_double)),
    ('temp', ctypes.POINTER(ctypes.c_double)),
    ('dens', ctypes.POINTER(ctypes.c_double)),
]

plasma_1D_data = struct_c__SA_plasma_1D_data
plasma_1D_init_offload = _libraries['libascot.so'].plasma_1D_init_offload
plasma_1D_init_offload.restype = ctypes.c_int32
plasma_1D_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1D_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
plasma_1D_free_offload = _libraries['libascot.so'].plasma_1D_free_offload
plasma_1D_free_offload.restype = None
plasma_1D_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1D_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
plasma_1D_init = _libraries['libascot.so'].plasma_1D_init
plasma_1D_init.restype = None
plasma_1D_init.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1D_data), ctypes.POINTER(struct_c__SA_plasma_1D_offload_data), ctypes.POINTER(ctypes.c_double)]
plasma_1D_eval_temp = _libraries['libascot.so'].plasma_1D_eval_temp
plasma_1D_eval_temp.restype = a5err
plasma_1D_eval_temp.argtypes = [ctypes.POINTER(ctypes.c_double), real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_plasma_1D_data)]
plasma_1D_eval_dens = _libraries['libascot.so'].plasma_1D_eval_dens
plasma_1D_eval_dens.restype = a5err
plasma_1D_eval_dens.argtypes = [ctypes.POINTER(ctypes.c_double), real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_plasma_1D_data)]
plasma_1D_eval_densandtemp = _libraries['libascot.so'].plasma_1D_eval_densandtemp
plasma_1D_eval_densandtemp.restype = a5err
plasma_1D_eval_densandtemp.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), real, ctypes.POINTER(struct_c__SA_plasma_1D_data)]
class struct_c__SA_plasma_1Dt_offload_data(Structure):
    pass

struct_c__SA_plasma_1Dt_offload_data._pack_ = 1 # source:False
struct_c__SA_plasma_1Dt_offload_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('n_time', ctypes.c_int32),
    ('n_species', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('mass', ctypes.c_double * 8),
    ('charge', ctypes.c_double * 8),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

plasma_1Dt_offload_data = struct_c__SA_plasma_1Dt_offload_data
class struct_c__SA_plasma_1Dt_data(Structure):
    pass

struct_c__SA_plasma_1Dt_data._pack_ = 1 # source:False
struct_c__SA_plasma_1Dt_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('n_time', ctypes.c_int32),
    ('n_species', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('mass', ctypes.c_double * 8),
    ('charge', ctypes.c_double * 8),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('rho', ctypes.POINTER(ctypes.c_double)),
    ('time', ctypes.POINTER(ctypes.c_double)),
    ('temp', ctypes.POINTER(ctypes.c_double)),
    ('dens', ctypes.POINTER(ctypes.c_double)),
]

plasma_1Dt_data = struct_c__SA_plasma_1Dt_data
plasma_1Dt_init_offload = _libraries['libascot.so'].plasma_1Dt_init_offload
plasma_1Dt_init_offload.restype = ctypes.c_int32
plasma_1Dt_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1Dt_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
plasma_1Dt_free_offload = _libraries['libascot.so'].plasma_1Dt_free_offload
plasma_1Dt_free_offload.restype = None
plasma_1Dt_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1Dt_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
plasma_1Dt_init = _libraries['libascot.so'].plasma_1Dt_init
plasma_1Dt_init.restype = None
plasma_1Dt_init.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1Dt_data), ctypes.POINTER(struct_c__SA_plasma_1Dt_offload_data), ctypes.POINTER(ctypes.c_double)]
plasma_1Dt_eval_temp = _libraries['libascot.so'].plasma_1Dt_eval_temp
plasma_1Dt_eval_temp.restype = a5err
plasma_1Dt_eval_temp.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_plasma_1Dt_data)]
plasma_1Dt_eval_dens = _libraries['libascot.so'].plasma_1Dt_eval_dens
plasma_1Dt_eval_dens.restype = a5err
plasma_1Dt_eval_dens.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_plasma_1Dt_data)]
plasma_1Dt_eval_densandtemp = _libraries['libascot.so'].plasma_1Dt_eval_densandtemp
plasma_1Dt_eval_densandtemp.restype = a5err
plasma_1Dt_eval_densandtemp.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), real, real, ctypes.POINTER(struct_c__SA_plasma_1Dt_data)]
class struct_c__SA_plasma_1DS_offload_data(Structure):
    pass

struct_c__SA_plasma_1DS_offload_data._pack_ = 1 # source:False
struct_c__SA_plasma_1DS_offload_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('rho_min', ctypes.c_double),
    ('rho_max', ctypes.c_double),
    ('n_species', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('mass', ctypes.c_double * 8),
    ('charge', ctypes.c_double * 8),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
]

plasma_1DS_offload_data = struct_c__SA_plasma_1DS_offload_data
class struct_c__SA_plasma_1DS_data(Structure):
    pass

struct_c__SA_plasma_1DS_data._pack_ = 1 # source:False
struct_c__SA_plasma_1DS_data._fields_ = [
    ('n_species', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('mass', ctypes.c_double * 8),
    ('charge', ctypes.c_double * 8),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('temp', struct_c__SA_interp1D_data * 2),
    ('dens', struct_c__SA_interp1D_data * 8),
]

plasma_1DS_data = struct_c__SA_plasma_1DS_data
plasma_1DS_init_offload = _libraries['libascot.so'].plasma_1DS_init_offload
plasma_1DS_init_offload.restype = ctypes.c_int32
plasma_1DS_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1DS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
plasma_1DS_free_offload = _libraries['libascot.so'].plasma_1DS_free_offload
plasma_1DS_free_offload.restype = None
plasma_1DS_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1DS_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
plasma_1DS_init = _libraries['libascot.so'].plasma_1DS_init
plasma_1DS_init.restype = None
plasma_1DS_init.argtypes = [ctypes.POINTER(struct_c__SA_plasma_1DS_data), ctypes.POINTER(struct_c__SA_plasma_1DS_offload_data), ctypes.POINTER(ctypes.c_double)]
plasma_1DS_eval_temp = _libraries['libascot.so'].plasma_1DS_eval_temp
plasma_1DS_eval_temp.restype = a5err
plasma_1DS_eval_temp.argtypes = [ctypes.POINTER(ctypes.c_double), real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_plasma_1DS_data)]
plasma_1DS_eval_dens = _libraries['libascot.so'].plasma_1DS_eval_dens
plasma_1DS_eval_dens.restype = a5err
plasma_1DS_eval_dens.argtypes = [ctypes.POINTER(ctypes.c_double), real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_plasma_1DS_data)]
plasma_1DS_eval_densandtemp = _libraries['libascot.so'].plasma_1DS_eval_densandtemp
plasma_1DS_eval_densandtemp.restype = a5err
plasma_1DS_eval_densandtemp.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), real, ctypes.POINTER(struct_c__SA_plasma_1DS_data)]

# values for enumeration 'plasma_type'
plasma_type__enumvalues = {
    0: 'plasma_type_1D',
    1: 'plasma_type_1Dt',
    2: 'plasma_type_1DS',
}
plasma_type_1D = 0
plasma_type_1Dt = 1
plasma_type_1DS = 2
plasma_type = ctypes.c_uint32 # enum
class struct_c__SA_plasma_offload_data(Structure):
    pass

struct_c__SA_plasma_offload_data._pack_ = 1 # source:False
struct_c__SA_plasma_offload_data._fields_ = [
    ('type', plasma_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('plasma_1D', plasma_1D_offload_data),
    ('plasma_1Dt', plasma_1Dt_offload_data),
    ('plasma_1DS', plasma_1DS_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

plasma_offload_data = struct_c__SA_plasma_offload_data
class struct_c__SA_plasma_data(Structure):
    pass

struct_c__SA_plasma_data._pack_ = 1 # source:False
struct_c__SA_plasma_data._fields_ = [
    ('type', plasma_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('plasma_1D', plasma_1D_data),
    ('plasma_1Dt', plasma_1Dt_data),
    ('plasma_1DS', plasma_1DS_data),
]

plasma_data = struct_c__SA_plasma_data
plasma_init_offload = _libraries['libascot.so'].plasma_init_offload
plasma_init_offload.restype = ctypes.c_int32
plasma_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_plasma_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
plasma_free_offload = _libraries['libascot.so'].plasma_free_offload
plasma_free_offload.restype = None
plasma_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_plasma_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
plasma_init = _libraries['libascot.so'].plasma_init
plasma_init.restype = ctypes.c_int32
plasma_init.argtypes = [ctypes.POINTER(struct_c__SA_plasma_data), ctypes.POINTER(struct_c__SA_plasma_offload_data), ctypes.POINTER(ctypes.c_double)]
plasma_eval_temp = _libraries['libascot.so'].plasma_eval_temp
plasma_eval_temp.restype = a5err
plasma_eval_temp.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, real, real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_plasma_data)]
plasma_eval_dens = _libraries['libascot.so'].plasma_eval_dens
plasma_eval_dens.restype = a5err
plasma_eval_dens.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, real, real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_plasma_data)]
plasma_eval_densandtemp = _libraries['libascot.so'].plasma_eval_densandtemp
plasma_eval_densandtemp.restype = a5err
plasma_eval_densandtemp.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), real, real, real, real, real, ctypes.POINTER(struct_c__SA_plasma_data)]
plasma_get_n_species = _libraries['libascot.so'].plasma_get_n_species
plasma_get_n_species.restype = ctypes.c_int32
plasma_get_n_species.argtypes = [ctypes.POINTER(struct_c__SA_plasma_data)]
plasma_get_species_mass = _libraries['libascot.so'].plasma_get_species_mass
plasma_get_species_mass.restype = ctypes.POINTER(ctypes.c_double)
plasma_get_species_mass.argtypes = [ctypes.POINTER(struct_c__SA_plasma_data)]
plasma_get_species_charge = _libraries['libascot.so'].plasma_get_species_charge
plasma_get_species_charge.restype = ctypes.POINTER(ctypes.c_double)
plasma_get_species_charge.argtypes = [ctypes.POINTER(struct_c__SA_plasma_data)]
plasma_get_species_znum = _libraries['libascot.so'].plasma_get_species_znum
plasma_get_species_znum.restype = ctypes.POINTER(ctypes.c_int32)
plasma_get_species_znum.argtypes = [ctypes.POINTER(struct_c__SA_plasma_data)]
plasma_get_species_anum = _libraries['libascot.so'].plasma_get_species_anum
plasma_get_species_anum.restype = ctypes.POINTER(ctypes.c_int32)
plasma_get_species_anum.argtypes = [ctypes.POINTER(struct_c__SA_plasma_data)]
class struct_c__SA_N0_1D_offload_data(Structure):
    pass

struct_c__SA_N0_1D_offload_data._pack_ = 1 # source:False
struct_c__SA_N0_1D_offload_data._fields_ = [
    ('n_rho', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('rho_min', ctypes.c_double),
    ('rho_max', ctypes.c_double),
    ('n_species', ctypes.c_int32),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('maxwellian', ctypes.c_int32 * 8),
    ('offload_array_length', ctypes.c_int32),
]

N0_1D_offload_data = struct_c__SA_N0_1D_offload_data
class struct_c__SA_N0_1D_data(Structure):
    pass

struct_c__SA_N0_1D_data._pack_ = 1 # source:False
struct_c__SA_N0_1D_data._fields_ = [
    ('n_species', ctypes.c_int32),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('maxwellian', ctypes.c_int32 * 8),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('n0', struct_c__SA_linint1D_data * 8),
    ('t0', struct_c__SA_linint1D_data * 8),
]

N0_1D_data = struct_c__SA_N0_1D_data
N0_1D_init_offload = _libraries['libascot.so'].N0_1D_init_offload
N0_1D_init_offload.restype = ctypes.c_int32
N0_1D_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_N0_1D_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
N0_1D_free_offload = _libraries['libascot.so'].N0_1D_free_offload
N0_1D_free_offload.restype = None
N0_1D_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_N0_1D_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
N0_1D_init = _libraries['libascot.so'].N0_1D_init
N0_1D_init.restype = None
N0_1D_init.argtypes = [ctypes.POINTER(struct_c__SA_N0_1D_data), ctypes.POINTER(struct_c__SA_N0_1D_offload_data), ctypes.POINTER(ctypes.c_double)]
N0_1D_eval_n0 = _libraries['libascot.so'].N0_1D_eval_n0
N0_1D_eval_n0.restype = a5err
N0_1D_eval_n0.argtypes = [ctypes.POINTER(ctypes.c_double), real, ctypes.POINTER(struct_c__SA_N0_1D_data)]
N0_1D_eval_t0 = _libraries['libascot.so'].N0_1D_eval_t0
N0_1D_eval_t0.restype = a5err
N0_1D_eval_t0.argtypes = [ctypes.POINTER(ctypes.c_double), real, ctypes.POINTER(struct_c__SA_N0_1D_data)]
N0_1D_get_n_species = _libraries['libascot.so'].N0_1D_get_n_species
N0_1D_get_n_species.restype = ctypes.c_int32
N0_1D_get_n_species.argtypes = [ctypes.POINTER(struct_c__SA_N0_1D_data)]
class struct_c__SA_N0_3D_offload_data(Structure):
    pass

struct_c__SA_N0_3D_offload_data._pack_ = 1 # source:False
struct_c__SA_N0_3D_offload_data._fields_ = [
    ('n_r', ctypes.c_int32),
    ('n_z', ctypes.c_int32),
    ('n_phi', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('r_min', ctypes.c_double),
    ('r_max', ctypes.c_double),
    ('z_min', ctypes.c_double),
    ('z_max', ctypes.c_double),
    ('phi_min', ctypes.c_double),
    ('phi_max', ctypes.c_double),
    ('n_species', ctypes.c_int32),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('maxwellian', ctypes.c_int32 * 8),
    ('offload_array_length', ctypes.c_int32),
]

N0_3D_offload_data = struct_c__SA_N0_3D_offload_data
class struct_c__SA_N0_3D_data(Structure):
    pass

class struct_c__SA_linint3D_data(Structure):
    pass

struct_c__SA_linint3D_data._pack_ = 1 # source:False
struct_c__SA_linint3D_data._fields_ = [
    ('n_x', ctypes.c_int32),
    ('n_y', ctypes.c_int32),
    ('n_z', ctypes.c_int32),
    ('bc_x', ctypes.c_int32),
    ('bc_y', ctypes.c_int32),
    ('bc_z', ctypes.c_int32),
    ('x_min', ctypes.c_double),
    ('x_max', ctypes.c_double),
    ('x_grid', ctypes.c_double),
    ('y_min', ctypes.c_double),
    ('y_max', ctypes.c_double),
    ('y_grid', ctypes.c_double),
    ('z_min', ctypes.c_double),
    ('z_max', ctypes.c_double),
    ('z_grid', ctypes.c_double),
    ('c', ctypes.POINTER(ctypes.c_double)),
]

struct_c__SA_N0_3D_data._pack_ = 1 # source:False
struct_c__SA_N0_3D_data._fields_ = [
    ('n_species', ctypes.c_int32),
    ('anum', ctypes.c_int32 * 8),
    ('znum', ctypes.c_int32 * 8),
    ('maxwellian', ctypes.c_int32 * 8),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('n0', struct_c__SA_linint3D_data * 8),
    ('t0', struct_c__SA_linint3D_data * 8),
]

N0_3D_data = struct_c__SA_N0_3D_data
N0_3D_init_offload = _libraries['libascot.so'].N0_3D_init_offload
N0_3D_init_offload.restype = ctypes.c_int32
N0_3D_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_N0_3D_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
N0_3D_free_offload = _libraries['libascot.so'].N0_3D_free_offload
N0_3D_free_offload.restype = None
N0_3D_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_N0_3D_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
N0_3D_init = _libraries['libascot.so'].N0_3D_init
N0_3D_init.restype = None
N0_3D_init.argtypes = [ctypes.POINTER(struct_c__SA_N0_3D_data), ctypes.POINTER(struct_c__SA_N0_3D_offload_data), ctypes.POINTER(ctypes.c_double)]
N0_3D_eval_n0 = _libraries['libascot.so'].N0_3D_eval_n0
N0_3D_eval_n0.restype = a5err
N0_3D_eval_n0.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, ctypes.POINTER(struct_c__SA_N0_3D_data)]
N0_3D_eval_t0 = _libraries['libascot.so'].N0_3D_eval_t0
N0_3D_eval_t0.restype = a5err
N0_3D_eval_t0.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, ctypes.POINTER(struct_c__SA_N0_3D_data)]
N0_3D_get_n_species = _libraries['libascot.so'].N0_3D_get_n_species
N0_3D_get_n_species.restype = ctypes.c_int32
N0_3D_get_n_species.argtypes = [ctypes.POINTER(struct_c__SA_N0_3D_data)]

# values for enumeration 'neutral_type'
neutral_type__enumvalues = {
    0: 'neutral_type_1D',
    1: 'neutral_type_3D',
}
neutral_type_1D = 0
neutral_type_3D = 1
neutral_type = ctypes.c_uint32 # enum
class struct_c__SA_neutral_offload_data(Structure):
    pass

struct_c__SA_neutral_offload_data._pack_ = 1 # source:False
struct_c__SA_neutral_offload_data._fields_ = [
    ('type', neutral_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('N01D', N0_1D_offload_data),
    ('N03D', N0_3D_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

neutral_offload_data = struct_c__SA_neutral_offload_data
class struct_c__SA_neutral_data(Structure):
    pass

struct_c__SA_neutral_data._pack_ = 1 # source:False
struct_c__SA_neutral_data._fields_ = [
    ('type', neutral_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('N01D', N0_1D_data),
    ('N03D', N0_3D_data),
]

neutral_data = struct_c__SA_neutral_data
neutral_init_offload = _libraries['libascot.so'].neutral_init_offload
neutral_init_offload.restype = ctypes.c_int32
neutral_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_neutral_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
neutral_free_offload = _libraries['libascot.so'].neutral_free_offload
neutral_free_offload.restype = None
neutral_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_neutral_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
neutral_init = _libraries['libascot.so'].neutral_init
neutral_init.restype = ctypes.c_int32
neutral_init.argtypes = [ctypes.POINTER(struct_c__SA_neutral_data), ctypes.POINTER(struct_c__SA_neutral_offload_data), ctypes.POINTER(ctypes.c_double)]
neutral_eval_n0 = _libraries['libascot.so'].neutral_eval_n0
neutral_eval_n0.restype = a5err
neutral_eval_n0.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, real, real, ctypes.POINTER(struct_c__SA_neutral_data)]
neutral_eval_t0 = _libraries['libascot.so'].neutral_eval_t0
neutral_eval_t0.restype = a5err
neutral_eval_t0.argtypes = [ctypes.POINTER(ctypes.c_double), real, real, real, real, real, ctypes.POINTER(struct_c__SA_neutral_data)]
neutral_get_n_species = _libraries['libascot.so'].neutral_get_n_species
neutral_get_n_species.restype = ctypes.c_int32
neutral_get_n_species.argtypes = [ctypes.POINTER(struct_c__SA_neutral_data)]
class struct_c__SA_wall_2d_offload_data(Structure):
    pass

struct_c__SA_wall_2d_offload_data._pack_ = 1 # source:False
struct_c__SA_wall_2d_offload_data._fields_ = [
    ('n', ctypes.c_int32),
    ('offload_array_length', ctypes.c_int32),
]

wall_2d_offload_data = struct_c__SA_wall_2d_offload_data
class struct_c__SA_wall_2d_data(Structure):
    pass

struct_c__SA_wall_2d_data._pack_ = 1 # source:False
struct_c__SA_wall_2d_data._fields_ = [
    ('n', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('wall_r', ctypes.POINTER(ctypes.c_double)),
    ('wall_z', ctypes.POINTER(ctypes.c_double)),
]

wall_2d_data = struct_c__SA_wall_2d_data
wall_2d_init_offload = _libraries['libascot.so'].wall_2d_init_offload
wall_2d_init_offload.restype = ctypes.c_int32
wall_2d_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_wall_2d_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
wall_2d_free_offload = _libraries['libascot.so'].wall_2d_free_offload
wall_2d_free_offload.restype = None
wall_2d_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_wall_2d_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
wall_2d_init = _libraries['libascot.so'].wall_2d_init
wall_2d_init.restype = None
wall_2d_init.argtypes = [ctypes.POINTER(struct_c__SA_wall_2d_data), ctypes.POINTER(struct_c__SA_wall_2d_offload_data), ctypes.POINTER(ctypes.c_double)]
wall_2d_inside = _libraries['libascot.so'].wall_2d_inside
wall_2d_inside.restype = ctypes.c_int32
wall_2d_inside.argtypes = [real, real, ctypes.POINTER(struct_c__SA_wall_2d_data)]
wall_2d_hit_wall = _libraries['libascot.so'].wall_2d_hit_wall
wall_2d_hit_wall.restype = ctypes.c_int32
wall_2d_hit_wall.argtypes = [real, real, real, real, real, real, ctypes.POINTER(struct_c__SA_wall_2d_data), ctypes.POINTER(ctypes.c_double)]
wall_2d_find_intersection = _libraries['libascot.so'].wall_2d_find_intersection
wall_2d_find_intersection.restype = ctypes.c_int32
wall_2d_find_intersection.argtypes = [real, real, real, real, ctypes.POINTER(struct_c__SA_wall_2d_data), ctypes.POINTER(ctypes.c_double)]
class struct_c__SA_wall_3d_offload_data(Structure):
    pass

struct_c__SA_wall_3d_offload_data._pack_ = 1 # source:False
struct_c__SA_wall_3d_offload_data._fields_ = [
    ('n', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('xmin', ctypes.c_double),
    ('xmax', ctypes.c_double),
    ('xgrid', ctypes.c_double),
    ('ymin', ctypes.c_double),
    ('ymax', ctypes.c_double),
    ('ygrid', ctypes.c_double),
    ('zmin', ctypes.c_double),
    ('zmax', ctypes.c_double),
    ('zgrid', ctypes.c_double),
    ('depth', ctypes.c_int32),
    ('ngrid', ctypes.c_int32),
    ('offload_array_length', ctypes.c_int32),
    ('int_offload_array_length', ctypes.c_int32),
]

wall_3d_offload_data = struct_c__SA_wall_3d_offload_data
class struct_c__SA_wall_3d_data(Structure):
    pass

struct_c__SA_wall_3d_data._pack_ = 1 # source:False
struct_c__SA_wall_3d_data._fields_ = [
    ('n', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('xmin', ctypes.c_double),
    ('xmax', ctypes.c_double),
    ('xgrid', ctypes.c_double),
    ('ymin', ctypes.c_double),
    ('ymax', ctypes.c_double),
    ('ygrid', ctypes.c_double),
    ('zmin', ctypes.c_double),
    ('zmax', ctypes.c_double),
    ('zgrid', ctypes.c_double),
    ('depth', ctypes.c_int32),
    ('ngrid', ctypes.c_int32),
    ('wall_tris', ctypes.POINTER(ctypes.c_double)),
    ('tree_array_size', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('tree_array', ctypes.POINTER(ctypes.c_int32)),
]

wall_3d_data = struct_c__SA_wall_3d_data
wall_3d_init_offload = _libraries['libascot.so'].wall_3d_init_offload
wall_3d_init_offload.restype = ctypes.c_int32
wall_3d_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_wall_3d_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32))]
wall_3d_free_offload = _libraries['libascot.so'].wall_3d_free_offload
wall_3d_free_offload.restype = None
wall_3d_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_wall_3d_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32))]
wall_3d_init_octree = _libraries['libascot.so'].wall_3d_init_octree
wall_3d_init_octree.restype = None
wall_3d_init_octree.argtypes = [ctypes.POINTER(struct_c__SA_wall_3d_offload_data), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32))]
wall_3d_init = _libraries['libascot.so'].wall_3d_init
wall_3d_init.restype = None
wall_3d_init.argtypes = [ctypes.POINTER(struct_c__SA_wall_3d_data), ctypes.POINTER(struct_c__SA_wall_3d_offload_data), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int32)]
wall_3d_hit_wall = _libraries['libascot.so'].wall_3d_hit_wall
wall_3d_hit_wall.restype = ctypes.c_int32
wall_3d_hit_wall.argtypes = [real, real, real, real, real, real, ctypes.POINTER(struct_c__SA_wall_3d_data), ctypes.POINTER(ctypes.c_double)]
wall_3d_hit_wall_full = _libraries['libascot.so'].wall_3d_hit_wall_full
wall_3d_hit_wall_full.restype = ctypes.c_int32
wall_3d_hit_wall_full.argtypes = [real, real, real, real, real, real, ctypes.POINTER(struct_c__SA_wall_3d_data), ctypes.POINTER(ctypes.c_double)]
wall_3d_tri_collision = _libraries['libascot.so'].wall_3d_tri_collision
wall_3d_tri_collision.restype = ctypes.c_double
wall_3d_tri_collision.argtypes = [ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3]
wall_3d_init_tree = _libraries['libascot.so'].wall_3d_init_tree
wall_3d_init_tree.restype = None
wall_3d_init_tree.argtypes = [ctypes.POINTER(struct_c__SA_wall_3d_data), ctypes.POINTER(ctypes.c_double)]
wall_3d_tri_in_cube = _libraries['libascot.so'].wall_3d_tri_in_cube
wall_3d_tri_in_cube.restype = ctypes.c_int32
wall_3d_tri_in_cube.argtypes = [ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3]
wall_3d_quad_collision = _libraries['libascot.so'].wall_3d_quad_collision
wall_3d_quad_collision.restype = ctypes.c_int32
wall_3d_quad_collision.argtypes = [ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3, ctypes.c_double * 3]

# values for enumeration 'wall_type'
wall_type__enumvalues = {
    0: 'wall_type_2D',
    1: 'wall_type_3D',
}
wall_type_2D = 0
wall_type_3D = 1
wall_type = ctypes.c_uint32 # enum
class struct_c__SA_wall_offload_data(Structure):
    pass

struct_c__SA_wall_offload_data._pack_ = 1 # source:False
struct_c__SA_wall_offload_data._fields_ = [
    ('type', wall_type),
    ('w2d', wall_2d_offload_data),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('w3d', wall_3d_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('int_offload_array_length', ctypes.c_int32),
]

wall_offload_data = struct_c__SA_wall_offload_data
class struct_c__SA_wall_data(Structure):
    pass

struct_c__SA_wall_data._pack_ = 1 # source:False
struct_c__SA_wall_data._fields_ = [
    ('type', wall_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('w2d', wall_2d_data),
    ('w3d', wall_3d_data),
]

wall_data = struct_c__SA_wall_data
wall_init_offload = _libraries['libascot.so'].wall_init_offload
wall_init_offload.restype = ctypes.c_int32
wall_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_wall_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32))]
wall_free_offload = _libraries['libascot.so'].wall_free_offload
wall_free_offload.restype = None
wall_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_wall_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32))]
wall_init = _libraries['libascot.so'].wall_init
wall_init.restype = ctypes.c_int32
wall_init.argtypes = [ctypes.POINTER(struct_c__SA_wall_data), ctypes.POINTER(struct_c__SA_wall_offload_data), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int32)]
wall_hit_wall = _libraries['libascot.so'].wall_hit_wall
wall_hit_wall.restype = ctypes.c_int32
wall_hit_wall.argtypes = [real, real, real, real, real, real, ctypes.POINTER(struct_c__SA_wall_data), ctypes.POINTER(ctypes.c_double)]
wall_get_n_elements = _libraries['libascot.so'].wall_get_n_elements
wall_get_n_elements.restype = ctypes.c_int32
wall_get_n_elements.argtypes = [ctypes.POINTER(struct_c__SA_wall_data)]
class struct_c__SA_boozer_offload_data(Structure):
    pass

struct_c__SA_boozer_offload_data._pack_ = 1 # source:False
struct_c__SA_boozer_offload_data._fields_ = [
    ('npsi', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('psi_min', ctypes.c_double),
    ('psi_max', ctypes.c_double),
    ('ntheta', ctypes.c_int32),
    ('nthetag', ctypes.c_int32),
    ('nrzs', ctypes.c_int32),
    ('offload_array_length', ctypes.c_int32),
]

boozer_offload_data = struct_c__SA_boozer_offload_data
class struct_c__SA_boozer_data(Structure):
    pass

struct_c__SA_boozer_data._pack_ = 1 # source:False
struct_c__SA_boozer_data._fields_ = [
    ('psi_min', ctypes.c_double),
    ('psi_max', ctypes.c_double),
    ('rs', ctypes.POINTER(ctypes.c_double)),
    ('zs', ctypes.POINTER(ctypes.c_double)),
    ('nrzs', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('nu_psitheta', struct_c__SA_interp2D_data),
    ('theta_psithetageom', struct_c__SA_interp2D_data),
]

boozer_data = struct_c__SA_boozer_data
boozer_init_offload = _libraries['libascot.so'].boozer_init_offload
boozer_init_offload.restype = ctypes.c_int32
boozer_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_boozer_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
boozer_free_offload = _libraries['libascot.so'].boozer_free_offload
boozer_free_offload.restype = None
boozer_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_boozer_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
boozer_init = _libraries['libascot.so'].boozer_init
boozer_init.restype = None
boozer_init.argtypes = [ctypes.POINTER(struct_c__SA_boozer_data), ctypes.POINTER(struct_c__SA_boozer_offload_data), ctypes.POINTER(ctypes.c_double)]
boozer_eval_psithetazeta = _libraries['libascot.so'].boozer_eval_psithetazeta
boozer_eval_psithetazeta.restype = a5err
boozer_eval_psithetazeta.argtypes = [ctypes.c_double * 12, ctypes.POINTER(ctypes.c_int32), real, real, real, ctypes.POINTER(struct_c__SA_B_field_data), ctypes.POINTER(struct_c__SA_boozer_data)]
class struct_c__SA_mhd_stat_offload_data(Structure):
    pass

struct_c__SA_mhd_stat_offload_data._pack_ = 1 # source:False
struct_c__SA_mhd_stat_offload_data._fields_ = [
    ('n_modes', ctypes.c_int32),
    ('nrho', ctypes.c_int32),
    ('rho_min', ctypes.c_double),
    ('rho_max', ctypes.c_double),
    ('nmode', ctypes.c_int32 * 512),
    ('mmode', ctypes.c_int32 * 512),
    ('amplitude_nm', ctypes.c_double * 512),
    ('omega_nm', ctypes.c_double * 512),
    ('phase_nm', ctypes.c_double * 512),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
]

mhd_stat_offload_data = struct_c__SA_mhd_stat_offload_data
class struct_c__SA_mhd_stat_data(Structure):
    pass

struct_c__SA_mhd_stat_data._pack_ = 1 # source:False
struct_c__SA_mhd_stat_data._fields_ = [
    ('n_modes', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('rho_min', ctypes.c_double),
    ('rho_max', ctypes.c_double),
    ('nmode', ctypes.c_int32 * 512),
    ('mmode', ctypes.c_int32 * 512),
    ('amplitude_nm', ctypes.c_double * 512),
    ('omega_nm', ctypes.c_double * 512),
    ('phase_nm', ctypes.c_double * 512),
    ('alpha_nm', struct_c__SA_interp1D_data * 512),
    ('phi_nm', struct_c__SA_interp1D_data * 512),
]

mhd_stat_data = struct_c__SA_mhd_stat_data
mhd_stat_init_offload = _libraries['libascot.so'].mhd_stat_init_offload
mhd_stat_init_offload.restype = ctypes.c_int32
mhd_stat_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_mhd_stat_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
mhd_stat_free_offload = _libraries['libascot.so'].mhd_stat_free_offload
mhd_stat_free_offload.restype = None
mhd_stat_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_mhd_stat_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
mhd_stat_init = _libraries['libascot.so'].mhd_stat_init
mhd_stat_init.restype = None
mhd_stat_init.argtypes = [ctypes.POINTER(struct_c__SA_mhd_stat_data), ctypes.POINTER(struct_c__SA_mhd_stat_offload_data), ctypes.POINTER(ctypes.c_double)]
mhd_stat_eval = _libraries['libascot.so'].mhd_stat_eval
mhd_stat_eval.restype = a5err
mhd_stat_eval.argtypes = [ctypes.c_double * 10, real, real, real, real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_boozer_data), ctypes.POINTER(struct_c__SA_mhd_stat_data), ctypes.POINTER(struct_c__SA_B_field_data)]
mhd_stat_perturbations = _libraries['libascot.so'].mhd_stat_perturbations
mhd_stat_perturbations.restype = a5err
mhd_stat_perturbations.argtypes = [ctypes.c_double * 7, real, real, real, real, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_boozer_data), ctypes.POINTER(struct_c__SA_mhd_stat_data), ctypes.POINTER(struct_c__SA_B_field_data)]
class struct_c__SA_mhd_nonstat_offload_data(Structure):
    pass

struct_c__SA_mhd_nonstat_offload_data._pack_ = 1 # source:False
struct_c__SA_mhd_nonstat_offload_data._fields_ = [
    ('n_modes', ctypes.c_int32),
    ('nrho', ctypes.c_int32),
    ('rho_min', ctypes.c_double),
    ('rho_max', ctypes.c_double),
    ('ntime', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('t_min', ctypes.c_double),
    ('t_max', ctypes.c_double),
    ('nmode', ctypes.c_int32 * 512),
    ('mmode', ctypes.c_int32 * 512),
    ('amplitude_nm', ctypes.c_double * 512),
    ('omega_nm', ctypes.c_double * 512),
    ('phase_nm', ctypes.c_double * 512),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

mhd_nonstat_offload_data = struct_c__SA_mhd_nonstat_offload_data
class struct_c__SA_mhd_nonstat_data(Structure):
    pass

struct_c__SA_mhd_nonstat_data._pack_ = 1 # source:False
struct_c__SA_mhd_nonstat_data._fields_ = [
    ('n_modes', ctypes.c_int32),
    ('nmode', ctypes.c_int32 * 512),
    ('mmode', ctypes.c_int32 * 512),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('amplitude_nm', ctypes.c_double * 512),
    ('omega_nm', ctypes.c_double * 512),
    ('phase_nm', ctypes.c_double * 512),
    ('alpha_nm', struct_c__SA_interp2D_data * 512),
    ('phi_nm', struct_c__SA_interp2D_data * 512),
]

mhd_nonstat_data = struct_c__SA_mhd_nonstat_data
mhd_nonstat_init_offload = _libraries['libascot.so'].mhd_nonstat_init_offload
mhd_nonstat_init_offload.restype = ctypes.c_int32
mhd_nonstat_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_mhd_nonstat_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
mhd_nonstat_free_offload = _libraries['libascot.so'].mhd_nonstat_free_offload
mhd_nonstat_free_offload.restype = None
mhd_nonstat_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_mhd_nonstat_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
mhd_nonstat_init = _libraries['libascot.so'].mhd_nonstat_init
mhd_nonstat_init.restype = None
mhd_nonstat_init.argtypes = [ctypes.POINTER(struct_c__SA_mhd_nonstat_data), ctypes.POINTER(struct_c__SA_mhd_nonstat_offload_data), ctypes.POINTER(ctypes.c_double)]
mhd_nonstat_eval = _libraries['libascot.so'].mhd_nonstat_eval
mhd_nonstat_eval.restype = a5err
mhd_nonstat_eval.argtypes = [ctypes.c_double * 10, real, real, real, real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_boozer_data), ctypes.POINTER(struct_c__SA_mhd_nonstat_data), ctypes.POINTER(struct_c__SA_B_field_data)]
mhd_nonstat_perturbations = _libraries['libascot.so'].mhd_nonstat_perturbations
mhd_nonstat_perturbations.restype = a5err
mhd_nonstat_perturbations.argtypes = [ctypes.c_double * 7, real, real, real, real, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_boozer_data), ctypes.POINTER(struct_c__SA_mhd_nonstat_data), ctypes.POINTER(struct_c__SA_B_field_data)]

# values for enumeration 'mhd_type'
mhd_type__enumvalues = {
    0: 'mhd_type_stat',
    1: 'mhd_type_nonstat',
}
mhd_type_stat = 0
mhd_type_nonstat = 1
mhd_type = ctypes.c_uint32 # enum
class struct_c__SA_mhd_offload_data(Structure):
    pass

struct_c__SA_mhd_offload_data._pack_ = 1 # source:False
struct_c__SA_mhd_offload_data._fields_ = [
    ('type', mhd_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('stat', mhd_stat_offload_data),
    ('nonstat', mhd_nonstat_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

mhd_offload_data = struct_c__SA_mhd_offload_data
class struct_c__SA_mhd_data(Structure):
    pass

struct_c__SA_mhd_data._pack_ = 1 # source:False
struct_c__SA_mhd_data._fields_ = [
    ('type', mhd_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('stat', mhd_stat_data),
    ('nonstat', mhd_nonstat_data),
]

mhd_data = struct_c__SA_mhd_data
mhd_init_offload = _libraries['libascot.so'].mhd_init_offload
mhd_init_offload.restype = ctypes.c_int32
mhd_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_mhd_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
mhd_free_offload = _libraries['libascot.so'].mhd_free_offload
mhd_free_offload.restype = None
mhd_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_mhd_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
mhd_init = _libraries['libascot.so'].mhd_init
mhd_init.restype = ctypes.c_int32
mhd_init.argtypes = [ctypes.POINTER(struct_c__SA_mhd_data), ctypes.POINTER(struct_c__SA_mhd_offload_data), ctypes.POINTER(ctypes.c_double)]
mhd_eval = _libraries['libascot.so'].mhd_eval
mhd_eval.restype = a5err
mhd_eval.argtypes = [ctypes.c_double * 10, real, real, real, real, ctypes.c_int32, ctypes.POINTER(struct_c__SA_boozer_data), ctypes.POINTER(struct_c__SA_mhd_data), ctypes.POINTER(struct_c__SA_B_field_data)]
mhd_perturbations = _libraries['libascot.so'].mhd_perturbations
mhd_perturbations.restype = a5err
mhd_perturbations.argtypes = [ctypes.c_double * 7, real, real, real, real, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_boozer_data), ctypes.POINTER(struct_c__SA_mhd_data), ctypes.POINTER(struct_c__SA_B_field_data)]
mhd_get_n_modes = _libraries['libascot.so'].mhd_get_n_modes
mhd_get_n_modes.restype = ctypes.c_int32
mhd_get_n_modes.argtypes = [ctypes.POINTER(struct_c__SA_mhd_data)]
mhd_get_nmode = _libraries['libascot.so'].mhd_get_nmode
mhd_get_nmode.restype = ctypes.POINTER(ctypes.c_int32)
mhd_get_nmode.argtypes = [ctypes.POINTER(struct_c__SA_mhd_data)]
mhd_get_mmode = _libraries['libascot.so'].mhd_get_mmode
mhd_get_mmode.restype = ctypes.POINTER(ctypes.c_int32)
mhd_get_mmode.argtypes = [ctypes.POINTER(struct_c__SA_mhd_data)]
mhd_get_amplitude = _libraries['libascot.so'].mhd_get_amplitude
mhd_get_amplitude.restype = ctypes.POINTER(ctypes.c_double)
mhd_get_amplitude.argtypes = [ctypes.POINTER(struct_c__SA_mhd_data)]
mhd_get_frequency = _libraries['libascot.so'].mhd_get_frequency
mhd_get_frequency.restype = ctypes.POINTER(ctypes.c_double)
mhd_get_frequency.argtypes = [ctypes.POINTER(struct_c__SA_mhd_data)]
mhd_get_phase = _libraries['libascot.so'].mhd_get_phase
mhd_get_phase.restype = ctypes.POINTER(ctypes.c_double)
mhd_get_phase.argtypes = [ctypes.POINTER(struct_c__SA_mhd_data)]
class struct_c__SA_asigma_loc_offload_data(Structure):
    pass

struct_c__SA_asigma_loc_offload_data._pack_ = 1 # source:False
struct_c__SA_asigma_loc_offload_data._fields_ = [
    ('N_reac', ctypes.c_int32),
    ('z_1', ctypes.c_int32 * 32),
    ('a_1', ctypes.c_int32 * 32),
    ('z_2', ctypes.c_int32 * 32),
    ('a_2', ctypes.c_int32 * 32),
    ('N_E', ctypes.c_int32 * 32),
    ('N_n', ctypes.c_int32 * 32),
    ('N_T', ctypes.c_int32 * 32),
    ('reac_type', ctypes.c_int32 * 32),
    ('offload_array_length', ctypes.c_int32),
]

asigma_loc_offload_data = struct_c__SA_asigma_loc_offload_data
class struct_c__SA_asigma_loc_data(Structure):
    pass

struct_c__SA_asigma_loc_data._pack_ = 1 # source:False
struct_c__SA_asigma_loc_data._fields_ = [
    ('N_reac', ctypes.c_int32),
    ('z_1', ctypes.c_int32 * 32),
    ('a_1', ctypes.c_int32 * 32),
    ('z_2', ctypes.c_int32 * 32),
    ('a_2', ctypes.c_int32 * 32),
    ('reac_type', ctypes.c_int32 * 32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('sigma', struct_c__SA_interp1D_data * 32),
    ('sigmav', struct_c__SA_interp2D_data * 32),
    ('BMSsigmav', struct_c__SA_interp3D_data * 32),
]

asigma_loc_data = struct_c__SA_asigma_loc_data
asigma_loc_init_offload = _libraries['libascot.so'].asigma_loc_init_offload
asigma_loc_init_offload.restype = ctypes.c_int32
asigma_loc_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_asigma_loc_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
asigma_loc_free_offload = _libraries['libascot.so'].asigma_loc_free_offload
asigma_loc_free_offload.restype = None
asigma_loc_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_asigma_loc_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
asigma_loc_init = _libraries['libascot.so'].asigma_loc_init
asigma_loc_init.restype = None
asigma_loc_init.argtypes = [ctypes.POINTER(struct_c__SA_asigma_loc_data), ctypes.POINTER(struct_c__SA_asigma_loc_offload_data), ctypes.POINTER(ctypes.c_double)]
asigma_loc_eval_sigma = _libraries['libascot.so'].asigma_loc_eval_sigma
asigma_loc_eval_sigma.restype = a5err
asigma_loc_eval_sigma.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, real, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_asigma_loc_data)]
asigma_loc_eval_sigmav = _libraries['libascot.so'].asigma_loc_eval_sigmav
asigma_loc_eval_sigmav.restype = a5err
asigma_loc_eval_sigmav.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.c_int32, real, ctypes.c_int32, ctypes.c_int32, real, real, real, real, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_asigma_loc_data)]

# values for enumeration 'asigma_type'
asigma_type__enumvalues = {
    0: 'asigma_type_loc',
}
asigma_type_loc = 0
asigma_type = ctypes.c_uint32 # enum

# values for enumeration 'asigma_reac_type'
asigma_reac_type__enumvalues = {
    1: 'sigma_ioniz',
    2: 'sigma_recomb',
    3: 'sigma_CX',
    4: 'sigmav_ioniz',
    5: 'sigmav_recomb',
    6: 'sigmav_CX',
    7: 'sigmav_BMS',
    8: 'sigmaveff_ioniz',
    9: 'sigmaveff_recomb',
    10: 'sigmaveff_CX',
}
sigma_ioniz = 1
sigma_recomb = 2
sigma_CX = 3
sigmav_ioniz = 4
sigmav_recomb = 5
sigmav_CX = 6
sigmav_BMS = 7
sigmaveff_ioniz = 8
sigmaveff_recomb = 9
sigmaveff_CX = 10
asigma_reac_type = ctypes.c_uint32 # enum
class struct_c__SA_asigma_offload_data(Structure):
    pass

struct_c__SA_asigma_offload_data._pack_ = 1 # source:False
struct_c__SA_asigma_offload_data._fields_ = [
    ('type', asigma_type),
    ('asigma_loc', asigma_loc_offload_data),
    ('offload_array_length', ctypes.c_int32),
]

asigma_offload_data = struct_c__SA_asigma_offload_data
class struct_c__SA_asigma_data(Structure):
    pass

struct_c__SA_asigma_data._pack_ = 1 # source:False
struct_c__SA_asigma_data._fields_ = [
    ('type', asigma_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('asigma_loc', asigma_loc_data),
]

asigma_data = struct_c__SA_asigma_data
asigma_init_offload = _libraries['libascot.so'].asigma_init_offload
asigma_init_offload.restype = ctypes.c_int32
asigma_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_asigma_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
asigma_free_offload = _libraries['libascot.so'].asigma_free_offload
asigma_free_offload.restype = None
asigma_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_asigma_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
asigma_extrapolate = _libraries['libascot.so'].asigma_extrapolate
asigma_extrapolate.restype = None
asigma_extrapolate.argtypes = [ctypes.c_int32]
asigma_init = _libraries['libascot.so'].asigma_init
asigma_init.restype = ctypes.c_int32
asigma_init.argtypes = [ctypes.POINTER(struct_c__SA_asigma_data), ctypes.POINTER(struct_c__SA_asigma_offload_data), ctypes.POINTER(ctypes.c_double)]
asigma_eval_sigma = _libraries['libascot.so'].asigma_eval_sigma
asigma_eval_sigma.restype = a5err
asigma_eval_sigma.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, real, asigma_reac_type, ctypes.POINTER(struct_c__SA_asigma_data)]
asigma_eval_sigmav = _libraries['libascot.so'].asigma_eval_sigmav
asigma_eval_sigmav.restype = a5err
asigma_eval_sigmav.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.c_int32, real, ctypes.c_int32, ctypes.c_int32, real, real, real, real, asigma_reac_type, ctypes.POINTER(struct_c__SA_asigma_data)]
class struct_c__SA_nbi_injector(Structure):
    pass

struct_c__SA_nbi_injector._pack_ = 1 # source:False
struct_c__SA_nbi_injector._fields_ = [
    ('id', ctypes.c_int32),
    ('n_beamlet', ctypes.c_int32),
    ('beamlet_x', ctypes.POINTER(ctypes.c_double)),
    ('beamlet_y', ctypes.POINTER(ctypes.c_double)),
    ('beamlet_z', ctypes.POINTER(ctypes.c_double)),
    ('beamlet_dx', ctypes.POINTER(ctypes.c_double)),
    ('beamlet_dy', ctypes.POINTER(ctypes.c_double)),
    ('beamlet_dz', ctypes.POINTER(ctypes.c_double)),
    ('power', ctypes.c_double),
    ('energy', ctypes.c_double),
    ('efrac', ctypes.c_double * 3),
    ('div_h', ctypes.c_double),
    ('div_v', ctypes.c_double),
    ('div_halo_frac', ctypes.c_double),
    ('div_halo_h', ctypes.c_double),
    ('div_halo_v', ctypes.c_double),
    ('anum', ctypes.c_int32),
    ('znum', ctypes.c_int32),
    ('mass', ctypes.c_double),
]

nbi_injector = struct_c__SA_nbi_injector
class struct_c__SA_nbi_offload_data(Structure):
    pass

struct_c__SA_nbi_offload_data._pack_ = 1 # source:False
struct_c__SA_nbi_offload_data._fields_ = [
    ('ninj', ctypes.c_int32),
    ('id', ctypes.c_int32 * 16),
    ('n_beamlet', ctypes.c_int32 * 16),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('power', ctypes.c_double * 16),
    ('energy', ctypes.c_double * 16),
    ('efrac', ctypes.c_double * 48),
    ('div_h', ctypes.c_double * 16),
    ('div_v', ctypes.c_double * 16),
    ('div_halo_frac', ctypes.c_double * 16),
    ('div_halo_h', ctypes.c_double * 16),
    ('div_halo_v', ctypes.c_double * 16),
    ('anum', ctypes.c_int32 * 16),
    ('znum', ctypes.c_int32 * 16),
    ('mass', ctypes.c_double * 16),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

nbi_offload_data = struct_c__SA_nbi_offload_data
class struct_c__SA_nbi_data(Structure):
    pass

struct_c__SA_nbi_data._pack_ = 1 # source:False
struct_c__SA_nbi_data._fields_ = [
    ('ninj', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('inj', struct_c__SA_nbi_injector * 16),
]

nbi_data = struct_c__SA_nbi_data
nbi_init_offload = _libraries['libascot.so'].nbi_init_offload
nbi_init_offload.restype = ctypes.c_int32
nbi_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_nbi_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
nbi_init = _libraries['libascot.so'].nbi_init
nbi_init.restype = None
nbi_init.argtypes = [ctypes.POINTER(struct_c__SA_nbi_data), ctypes.POINTER(struct_c__SA_nbi_offload_data), ctypes.POINTER(ctypes.c_double)]
nbi_free_offload = _libraries['libascot.so'].nbi_free_offload
nbi_free_offload.restype = None
nbi_free_offload.argtypes = [ctypes.POINTER(struct_c__SA_nbi_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
nbi_inject = _libraries['libascot.so'].nbi_inject
nbi_inject.restype = None
nbi_inject.argtypes = [ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(struct_c__SA_nbi_injector), ctypes.POINTER(ctypes.POINTER(None))]

# values for enumeration 'SIMULATION_MODE'
SIMULATION_MODE__enumvalues = {
    1: 'simulate_mode_fo',
    2: 'simulate_mode_gc',
    3: 'simulate_mode_hybrid',
    4: 'simulate_mode_ml',
}
simulate_mode_fo = 1
simulate_mode_gc = 2
simulate_mode_hybrid = 3
simulate_mode_ml = 4
SIMULATION_MODE = ctypes.c_uint32 # enum
class struct_c__SA_sim_offload_data(Structure):
    pass

struct_c__SA_sim_offload_data._pack_ = 1 # source:False
struct_c__SA_sim_offload_data._fields_ = [
    ('B_offload_data', B_field_offload_data),
    ('E_offload_data', E_field_offload_data),
    ('plasma_offload_data', plasma_offload_data),
    ('neutral_offload_data', neutral_offload_data),
    ('wall_offload_data', wall_offload_data),
    ('boozer_offload_data', boozer_offload_data),
    ('mhd_offload_data', mhd_offload_data),
    ('asigma_offload_data', asigma_offload_data),
    ('nbi_offload_data', nbi_offload_data),
    ('diag_offload_data', diag_offload_data),
    ('sim_mode', ctypes.c_int32),
    ('enable_ada', ctypes.c_int32),
    ('record_mode', ctypes.c_int32),
    ('fix_usrdef_use', ctypes.c_int32),
    ('fix_usrdef_val', ctypes.c_double),
    ('fix_gyrodef_nstep', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('ada_tol_orbfol', ctypes.c_double),
    ('ada_tol_clmbcol', ctypes.c_double),
    ('ada_max_drho', ctypes.c_double),
    ('ada_max_dphi', ctypes.c_double),
    ('enable_orbfol', ctypes.c_int32),
    ('enable_clmbcol', ctypes.c_int32),
    ('enable_mhd', ctypes.c_int32),
    ('enable_atomic', ctypes.c_int32),
    ('disable_gctransform', ctypes.c_int32),
    ('disable_energyccoll', ctypes.c_int32),
    ('disable_pitchccoll', ctypes.c_int32),
    ('disable_gcdiffccoll', ctypes.c_int32),
    ('reverse_time', ctypes.c_int32),
    ('endcond_active', ctypes.c_int32),
    ('endcond_lim_simtime', ctypes.c_double),
    ('endcond_max_mileage', ctypes.c_double),
    ('endcond_max_cputime', ctypes.c_double),
    ('endcond_min_rho', ctypes.c_double),
    ('endcond_max_rho', ctypes.c_double),
    ('endcond_min_ekin', ctypes.c_double),
    ('endcond_min_thermal', ctypes.c_double),
    ('endcond_max_tororb', ctypes.c_double),
    ('endcond_max_polorb', ctypes.c_double),
    ('endcond_torandpol', ctypes.c_int32),
    ('hdf5_in', ctypes.c_char * 256),
    ('hdf5_out', ctypes.c_char * 256),
    ('qid', ctypes.c_char * 256),
    ('description', ctypes.c_char * 256),
    ('mpi_rank', ctypes.c_int32),
    ('mpi_size', ctypes.c_int32),
    ('qid_options', ctypes.c_char * 256),
    ('qid_bfield', ctypes.c_char * 256),
    ('qid_efield', ctypes.c_char * 256),
    ('qid_marker', ctypes.c_char * 256),
    ('qid_wall', ctypes.c_char * 256),
    ('qid_plasma', ctypes.c_char * 256),
    ('qid_neutral', ctypes.c_char * 256),
    ('qid_boozer', ctypes.c_char * 256),
    ('qid_mhd', ctypes.c_char * 256),
    ('qid_asigma', ctypes.c_char * 256),
    ('qid_nbi', ctypes.c_char * 256),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

sim_offload_data = struct_c__SA_sim_offload_data
class struct_c__SA_sim_data(Structure):
    pass

class struct_c__SA_mccc_data(Structure):
    pass

struct_c__SA_mccc_data._pack_ = 1 # source:False
struct_c__SA_mccc_data._fields_ = [
    ('usetabulated', ctypes.c_int32),
    ('include_energy', ctypes.c_int32),
    ('include_pitch', ctypes.c_int32),
    ('include_gcdiff', ctypes.c_int32),
]

struct_c__SA_sim_data._pack_ = 1 # source:False
struct_c__SA_sim_data._fields_ = [
    ('B_data', B_field_data),
    ('E_data', E_field_data),
    ('plasma_data', plasma_data),
    ('neutral_data', neutral_data),
    ('wall_data', wall_data),
    ('boozer_data', boozer_data),
    ('mhd_data', mhd_data),
    ('asigma_data', asigma_data),
    ('nbi_data', nbi_data),
    ('diag_data', diag_data),
    ('random_data', ctypes.POINTER(None)),
    ('mccc_data', struct_c__SA_mccc_data),
    ('sim_mode', ctypes.c_int32),
    ('enable_ada', ctypes.c_int32),
    ('record_mode', ctypes.c_int32),
    ('fix_usrdef_use', ctypes.c_int32),
    ('fix_usrdef_val', ctypes.c_double),
    ('fix_gyrodef_nstep', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('ada_tol_orbfol', ctypes.c_double),
    ('ada_tol_clmbcol', ctypes.c_double),
    ('ada_max_drho', ctypes.c_double),
    ('ada_max_dphi', ctypes.c_double),
    ('enable_orbfol', ctypes.c_int32),
    ('enable_clmbcol', ctypes.c_int32),
    ('enable_mhd', ctypes.c_int32),
    ('enable_atomic', ctypes.c_int32),
    ('disable_gctransform', ctypes.c_int32),
    ('disable_energyccoll', ctypes.c_int32),
    ('disable_pitchccoll', ctypes.c_int32),
    ('disable_gcdiffccoll', ctypes.c_int32),
    ('reverse_time', ctypes.c_int32),
    ('endcond_active', ctypes.c_int32),
    ('endcond_lim_simtime', ctypes.c_double),
    ('endcond_max_mileage', ctypes.c_double),
    ('endcond_max_cputime', ctypes.c_double),
    ('endcond_min_rho', ctypes.c_double),
    ('endcond_max_rho', ctypes.c_double),
    ('endcond_min_ekin', ctypes.c_double),
    ('endcond_min_thermal', ctypes.c_double),
    ('endcond_max_tororb', ctypes.c_double),
    ('endcond_max_polorb', ctypes.c_double),
    ('endcond_torandpol', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

sim_data = struct_c__SA_sim_data
simulate_init_offload = _libraries['libascot.so'].simulate_init_offload
simulate_init_offload.restype = None
simulate_init_offload.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data)]
sim_init = _libraries['libascot.so'].sim_init
sim_init.restype = None
sim_init.argtypes = [ctypes.POINTER(struct_c__SA_sim_data), ctypes.POINTER(struct_c__SA_sim_offload_data)]
class struct_c__SA_offload_package(Structure):
    pass

struct_c__SA_offload_package._pack_ = 1 # source:False
struct_c__SA_offload_package._fields_ = [
    ('offload_array_length', ctypes.c_uint64),
    ('unpack_pos', ctypes.c_uint64),
    ('int_offload_array_length', ctypes.c_uint64),
    ('int_unpack_pos', ctypes.c_uint64),
]

simulate = _libraries['libascot.so'].simulate
simulate.restype = None
simulate.argtypes = [ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.POINTER(struct_c__SA_offload_package), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_double)]
mpi_interface_init = _libraries['libascot.so'].mpi_interface_init
mpi_interface_init.restype = None
mpi_interface_init.argtypes = [ctypes.c_int32, ctypes.POINTER(ctypes.POINTER(ctypes.c_char)), ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32)]
mpi_interface_finalize = _libraries['libascot.so'].mpi_interface_finalize
mpi_interface_finalize.restype = None
mpi_interface_finalize.argtypes = []
mpi_my_particles = _libraries['libascot.so'].mpi_my_particles
mpi_my_particles.restype = None
mpi_my_particles.argtypes = [ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32]
mpi_gather_particlestate = _libraries['libascot.so'].mpi_gather_particlestate
mpi_gather_particlestate.restype = None
mpi_gather_particlestate.argtypes = [ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(ctypes.POINTER(struct_c__SA_particle_state)), ctypes.POINTER(ctypes.c_int32), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32]
mpi_gather_diag = _libraries['libascot.so'].mpi_gather_diag
mpi_gather_diag.restype = None
mpi_gather_diag.argtypes = [ctypes.POINTER(struct_c__SA_diag_offload_data), ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32]

# values for enumeration 'ENDCOND_FLAG'
ENDCOND_FLAG__enumvalues = {
    1: 'endcond_tlim',
    2: 'endcond_emin',
    4: 'endcond_therm',
    8: 'endcond_wall',
    16: 'endcond_rhomin',
    32: 'endcond_rhomax',
    64: 'endcond_polmax',
    128: 'endcond_tormax',
    256: 'endcond_cpumax',
    512: 'endcond_hybrid',
    1024: 'endcond_neutr',
    2048: 'endcond_ioniz',
}
endcond_tlim = 1
endcond_emin = 2
endcond_therm = 4
endcond_wall = 8
endcond_rhomin = 16
endcond_rhomax = 32
endcond_polmax = 64
endcond_tormax = 128
endcond_cpumax = 256
endcond_hybrid = 512
endcond_neutr = 1024
endcond_ioniz = 2048
ENDCOND_FLAG = ctypes.c_uint32 # enum
endcond_check_gc = _libraries['libascot.so'].endcond_check_gc
endcond_check_gc.restype = None
endcond_check_gc.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_particle_simd_gc), ctypes.POINTER(struct_c__SA_sim_data)]
endcond_check_fo = _libraries['libascot.so'].endcond_check_fo
endcond_check_fo.restype = None
endcond_check_fo.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_particle_simd_fo), ctypes.POINTER(struct_c__SA_sim_data)]
endcond_check_ml = _libraries['libascot.so'].endcond_check_ml
endcond_check_ml.restype = None
endcond_check_ml.argtypes = [ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.POINTER(struct_c__SA_particle_simd_ml), ctypes.POINTER(struct_c__SA_sim_data)]
endcond_parse = _libraries['libascot.so'].endcond_parse
endcond_parse.restype = None
endcond_parse.argtypes = [ctypes.c_int32, ctypes.POINTER(ctypes.c_int32)]
endcond_parse2str = _libraries['libascot.so'].endcond_parse2str
endcond_parse2str.restype = None
endcond_parse2str.argtypes = [ctypes.c_int32, ctypes.POINTER(ctypes.c_char)]

# values for enumeration 'input_group'
input_group__enumvalues = {
    1: 'hdf5_input_options',
    2: 'hdf5_input_bfield',
    4: 'hdf5_input_efield',
    8: 'hdf5_input_plasma',
    16: 'hdf5_input_neutral',
    32: 'hdf5_input_wall',
    64: 'hdf5_input_marker',
    128: 'hdf5_input_boozer',
    256: 'hdf5_input_mhd',
    512: 'hdf5_input_asigma',
    1024: 'hdf5_input_nbi',
}
hdf5_input_options = 1
hdf5_input_bfield = 2
hdf5_input_efield = 4
hdf5_input_plasma = 8
hdf5_input_neutral = 16
hdf5_input_wall = 32
hdf5_input_marker = 64
hdf5_input_boozer = 128
hdf5_input_mhd = 256
hdf5_input_asigma = 512
hdf5_input_nbi = 1024
input_group = ctypes.c_uint32 # enum
hdf5_interface_read_input = _libraries['libascot.so'].hdf5_interface_read_input
hdf5_interface_read_input.restype = ctypes.c_int32
hdf5_interface_read_input.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.c_int32, ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(struct_c__SA_input_particle)), ctypes.POINTER(ctypes.c_int32)]
hdf5_interface_init_results = _libraries['libascot.so'].hdf5_interface_init_results
hdf5_interface_init_results.restype = ctypes.c_int32
hdf5_interface_init_results.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_char)]
hdf5_interface_write_state = _libraries['libascot.so'].hdf5_interface_write_state
hdf5_interface_write_state.restype = ctypes.c_int32
hdf5_interface_write_state.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_char), integer, ctypes.POINTER(struct_c__SA_particle_state)]
hdf5_interface_write_diagnostics = _libraries['libascot.so'].hdf5_interface_write_diagnostics
hdf5_interface_write_diagnostics.restype = ctypes.c_int32
hdf5_interface_write_diagnostics.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_char)]
hdf5_generate_qid = _libraries['libascot.so'].hdf5_generate_qid
hdf5_generate_qid.restype = None
hdf5_generate_qid.argtypes = [ctypes.POINTER(ctypes.c_char)]
libascot_allocate_input_particles = _libraries['libascot.so'].libascot_allocate_input_particles
libascot_allocate_input_particles.restype = ctypes.POINTER(struct_c__SA_input_particle)
libascot_allocate_input_particles.argtypes = [ctypes.c_int32]
libascot_allocate_particle_states = _libraries['libascot.so'].libascot_allocate_particle_states
libascot_allocate_particle_states.restype = ctypes.POINTER(struct_c__SA_particle_state)
libascot_allocate_particle_states.argtypes = [ctypes.c_int32]
libascot_allocate_reals = _libraries['libascot.so'].libascot_allocate_reals
libascot_allocate_reals.restype = ctypes.POINTER(ctypes.c_double)
libascot_allocate_reals.argtypes = [size_t]
libascot_deallocate = _libraries['libascot.so'].libascot_deallocate
libascot_deallocate.restype = None
libascot_deallocate.argtypes = [ctypes.POINTER(None)]

# values for enumeration 'Reaction'
Reaction__enumvalues = {
    1: 'DT_He4n',
    2: 'DHe3_He4p',
    3: 'DD_Tp',
    4: 'DD_He3n',
}
DT_He4n = 1
DHe3_He4p = 2
DD_Tp = 3
DD_He3n = 4
Reaction = ctypes.c_uint32 # enum
boschhale_reaction = _libraries['libascot.so'].boschhale_reaction
boschhale_reaction.restype = None
boschhale_reaction.argtypes = [Reaction, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
boschhale_sigma = _libraries['libascot.so'].boschhale_sigma
boschhale_sigma.restype = real
boschhale_sigma.argtypes = [Reaction, real]
boschhale_sigmav = _libraries['libascot.so'].boschhale_sigmav
boschhale_sigmav.restype = real
boschhale_sigmav.argtypes = [Reaction, real]
class struct_c__SA_afsi_thermal_data(Structure):
    pass

struct_c__SA_afsi_thermal_data._pack_ = 1 # source:False
struct_c__SA_afsi_thermal_data._fields_ = [
    ('n_r', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('min_r', ctypes.c_double),
    ('max_r', ctypes.c_double),
    ('n_phi', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('min_phi', ctypes.c_double),
    ('max_phi', ctypes.c_double),
    ('n_z', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('min_z', ctypes.c_double),
    ('max_z', ctypes.c_double),
    ('temperature', ctypes.POINTER(ctypes.c_double)),
    ('density', ctypes.POINTER(ctypes.c_double)),
]

afsi_thermal_data = struct_c__SA_afsi_thermal_data
class struct_c__SA_afsi_data(Structure):
    pass

struct_c__SA_afsi_data._pack_ = 1 # source:False
struct_c__SA_afsi_data._fields_ = [
    ('type', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('dist_5D', ctypes.POINTER(struct_c__SA_dist_5D_data)),
    ('dist_thermal', ctypes.POINTER(struct_c__SA_afsi_thermal_data)),
]

afsi_data = struct_c__SA_afsi_data
afsi_run = _libraries['libascot.so'].afsi_run
afsi_run.restype = None
afsi_run.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), Reaction, ctypes.c_int32, ctypes.POINTER(struct_c__SA_afsi_data), ctypes.POINTER(struct_c__SA_afsi_data), real, ctypes.POINTER(struct_c__SA_dist_5D_offload_data), ctypes.POINTER(struct_c__SA_dist_5D_offload_data), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
afsi_test_dist = _libraries['libascot.so'].afsi_test_dist
afsi_test_dist.restype = None
afsi_test_dist.argtypes = [ctypes.POINTER(struct_c__SA_dist_5D_data)]
afsi_test_thermal = _libraries['libascot.so'].afsi_test_thermal
afsi_test_thermal.restype = None
afsi_test_thermal.argtypes = []
pack_offload_array = _libraries['libascot.so'].pack_offload_array
pack_offload_array.restype = ctypes.c_int32
pack_offload_array.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.POINTER(struct_c__SA_offload_package), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_int32))]
prepare_markers = _libraries['libascot.so'].prepare_markers
prepare_markers.restype = ctypes.c_int32
prepare_markers.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(ctypes.POINTER(struct_c__SA_input_particle)), ctypes.POINTER(ctypes.POINTER(struct_c__SA_particle_state)), ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_double)]
write_rungroup = _libraries['libascot.so'].write_rungroup
write_rungroup.restype = ctypes.c_int32
write_rungroup.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(ctypes.c_char)]
offload_and_simulate = _libraries['libascot.so'].offload_and_simulate
offload_and_simulate.restype = ctypes.c_int32
offload_and_simulate.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_state), ctypes.POINTER(struct_c__SA_offload_package), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.POINTER(struct_c__SA_particle_state)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
write_output = _libraries['libascot.so'].write_output
write_output.restype = ctypes.c_int32
write_output.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(struct_c__SA_particle_state), ctypes.c_int32, ctypes.POINTER(ctypes.c_double)]
print_marker_summary = _libraries['libascot.so'].print_marker_summary
print_marker_summary.restype = None
print_marker_summary.argtypes = [ctypes.POINTER(struct_c__SA_particle_state), ctypes.c_int32]
biosaw_calc_B = _libraries['libascot.so'].biosaw_calc_B
biosaw_calc_B.restype = None
biosaw_calc_B.argtypes = [ctypes.c_int32, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int32, ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
__all__ = \
    ['B_2DS_data', 'B_2DS_eval_B', 'B_2DS_eval_B_dB',
    'B_2DS_eval_psi', 'B_2DS_eval_psi_dpsi', 'B_2DS_eval_rho_drho',
    'B_2DS_free_offload', 'B_2DS_get_axis_rz', 'B_2DS_init',
    'B_2DS_init_offload', 'B_2DS_offload_data', 'B_3DS_data',
    'B_3DS_eval_B', 'B_3DS_eval_B_dB', 'B_3DS_eval_psi',
    'B_3DS_eval_psi_dpsi', 'B_3DS_eval_rho_drho',
    'B_3DS_free_offload', 'B_3DS_get_axis_rz', 'B_3DS_init',
    'B_3DS_init_offload', 'B_3DS_offload_data', 'B_GS_data',
    'B_GS_eval_B', 'B_GS_eval_B_dB', 'B_GS_eval_psi',
    'B_GS_eval_psi_dpsi', 'B_GS_eval_rho_drho', 'B_GS_free_offload',
    'B_GS_get_axis_rz', 'B_GS_init', 'B_GS_init_offload',
    'B_GS_offload_data', 'B_STS_data', 'B_STS_eval_B',
    'B_STS_eval_B_dB', 'B_STS_eval_psi', 'B_STS_eval_psi_dpsi',
    'B_STS_eval_rho_drho', 'B_STS_free_offload', 'B_STS_get_axis_rz',
    'B_STS_init', 'B_STS_init_offload', 'B_STS_offload_data',
    'B_TC_data', 'B_TC_eval_B', 'B_TC_eval_B_dB', 'B_TC_eval_psi',
    'B_TC_eval_psi_dpsi', 'B_TC_eval_rho_drho', 'B_TC_free_offload',
    'B_TC_get_axis_rz', 'B_TC_init', 'B_TC_init_offload',
    'B_TC_offload_data', 'B_field_data', 'B_field_eval_B',
    'B_field_eval_B_dB', 'B_field_eval_psi', 'B_field_eval_psi_dpsi',
    'B_field_eval_rho', 'B_field_eval_rho_drho',
    'B_field_free_offload', 'B_field_get_axis_rz', 'B_field_init',
    'B_field_init_offload', 'B_field_offload_data', 'B_field_type',
    'B_field_type_2DS', 'B_field_type_3DS', 'B_field_type_GS',
    'B_field_type_STS', 'B_field_type_TC', 'DD_He3n', 'DD_Tp',
    'DHe3_He4p', 'DT_He4n', 'ENDCOND_FLAG', 'E_1DS_data',
    'E_1DS_eval_E', 'E_1DS_free_offload', 'E_1DS_init',
    'E_1DS_init_offload', 'E_1DS_offload_data', 'E_TC_data',
    'E_TC_eval_E', 'E_TC_free_offload', 'E_TC_init',
    'E_TC_init_offload', 'E_TC_offload_data', 'E_field_data',
    'E_field_eval_E', 'E_field_free_offload', 'E_field_init',
    'E_field_init_offload', 'E_field_offload_data', 'E_field_type',
    'E_field_type_1DS', 'E_field_type_TC', 'N0_1D_data',
    'N0_1D_eval_n0', 'N0_1D_eval_t0', 'N0_1D_free_offload',
    'N0_1D_get_n_species', 'N0_1D_init', 'N0_1D_init_offload',
    'N0_1D_offload_data', 'N0_3D_data', 'N0_3D_eval_n0',
    'N0_3D_eval_t0', 'N0_3D_free_offload', 'N0_3D_get_n_species',
    'N0_3D_init', 'N0_3D_init_offload', 'N0_3D_offload_data',
    'Reaction', 'SIMULATION_MODE', 'a5err', 'afsi_data', 'afsi_run',
    'afsi_test_dist', 'afsi_test_thermal', 'afsi_thermal_data',
    'asigma_data', 'asigma_eval_sigma', 'asigma_eval_sigmav',
    'asigma_extrapolate', 'asigma_free_offload', 'asigma_init',
    'asigma_init_offload', 'asigma_loc_data', 'asigma_loc_eval_sigma',
    'asigma_loc_eval_sigmav', 'asigma_loc_free_offload',
    'asigma_loc_init', 'asigma_loc_init_offload',
    'asigma_loc_offload_data', 'asigma_offload_data',
    'asigma_reac_type', 'asigma_type', 'asigma_type_loc',
    'biosaw_calc_B', 'boozer_data', 'boozer_eval_psithetazeta',
    'boozer_free_offload', 'boozer_init', 'boozer_init_offload',
    'boozer_offload_data', 'boschhale_reaction', 'boschhale_sigma',
    'boschhale_sigmav', 'diag_data', 'diag_free', 'diag_free_offload',
    'diag_init', 'diag_init_offload', 'diag_offload_data',
    'diag_orb_check_plane_crossing', 'diag_orb_check_radial_crossing',
    'diag_orb_data', 'diag_orb_free', 'diag_orb_init',
    'diag_orb_offload_data', 'diag_orb_update_fo',
    'diag_orb_update_gc', 'diag_orb_update_ml', 'diag_sum',
    'diag_transcoef_data', 'diag_transcoef_free',
    'diag_transcoef_init', 'diag_transcoef_link',
    'diag_transcoef_offload_data', 'diag_transcoef_update_fo',
    'diag_transcoef_update_gc', 'diag_transcoef_update_ml',
    'diag_update_fo', 'diag_update_gc', 'diag_update_ml',
    'dist_5D_data', 'dist_5D_index', 'dist_5D_init',
    'dist_5D_offload_data', 'dist_5D_update_fo', 'dist_5D_update_gc',
    'dist_6D_data', 'dist_6D_init', 'dist_6D_offload_data',
    'dist_6D_update_fo', 'dist_6D_update_gc', 'dist_COM_data',
    'dist_COM_init', 'dist_COM_offload_data', 'dist_COM_update_fo',
    'dist_COM_update_gc', 'dist_rho5D_data', 'dist_rho5D_init',
    'dist_rho5D_offload_data', 'dist_rho5D_update_fo',
    'dist_rho5D_update_gc', 'dist_rho6D_data', 'dist_rho6D_init',
    'dist_rho6D_offload_data', 'dist_rho6D_update_fo',
    'dist_rho6D_update_gc', 'endcond_check_fo', 'endcond_check_gc',
    'endcond_check_ml', 'endcond_cpumax', 'endcond_emin',
    'endcond_hybrid', 'endcond_ioniz', 'endcond_neutr',
    'endcond_parse', 'endcond_parse2str', 'endcond_polmax',
    'endcond_rhomax', 'endcond_rhomin', 'endcond_therm',
    'endcond_tlim', 'endcond_tormax', 'endcond_wall',
    'hdf5_generate_qid', 'hdf5_input_asigma', 'hdf5_input_bfield',
    'hdf5_input_boozer', 'hdf5_input_efield', 'hdf5_input_marker',
    'hdf5_input_mhd', 'hdf5_input_nbi', 'hdf5_input_neutral',
    'hdf5_input_options', 'hdf5_input_plasma', 'hdf5_input_wall',
    'hdf5_interface_init_results', 'hdf5_interface_read_input',
    'hdf5_interface_write_diagnostics', 'hdf5_interface_write_state',
    'input_group', 'input_particle', 'input_particle_type',
    'input_particle_type_gc', 'input_particle_type_ml',
    'input_particle_type_p', 'input_particle_type_s', 'integer',
    'libascot_allocate_input_particles',
    'libascot_allocate_particle_states', 'libascot_allocate_reals',
    'libascot_deallocate', 'mhd_data', 'mhd_eval', 'mhd_free_offload',
    'mhd_get_amplitude', 'mhd_get_frequency', 'mhd_get_mmode',
    'mhd_get_n_modes', 'mhd_get_nmode', 'mhd_get_phase', 'mhd_init',
    'mhd_init_offload', 'mhd_nonstat_data', 'mhd_nonstat_eval',
    'mhd_nonstat_free_offload', 'mhd_nonstat_init',
    'mhd_nonstat_init_offload', 'mhd_nonstat_offload_data',
    'mhd_nonstat_perturbations', 'mhd_offload_data',
    'mhd_perturbations', 'mhd_stat_data', 'mhd_stat_eval',
    'mhd_stat_free_offload', 'mhd_stat_init', 'mhd_stat_init_offload',
    'mhd_stat_offload_data', 'mhd_stat_perturbations', 'mhd_type',
    'mhd_type_nonstat', 'mhd_type_stat', 'mpi_gather_diag',
    'mpi_gather_particlestate', 'mpi_interface_finalize',
    'mpi_interface_init', 'mpi_my_particles', 'nbi_data',
    'nbi_free_offload', 'nbi_init', 'nbi_init_offload', 'nbi_inject',
    'nbi_injector', 'nbi_offload_data', 'neutral_data',
    'neutral_eval_n0', 'neutral_eval_t0', 'neutral_free_offload',
    'neutral_get_n_species', 'neutral_init', 'neutral_init_offload',
    'neutral_offload_data', 'neutral_type', 'neutral_type_1D',
    'neutral_type_3D', 'offload_and_simulate', 'pack_offload_array',
    'particle', 'particle_copy_fo', 'particle_copy_gc',
    'particle_copy_ml', 'particle_cycle_fo', 'particle_cycle_gc',
    'particle_cycle_ml', 'particle_fo_to_gc', 'particle_fo_to_state',
    'particle_gc', 'particle_gc_to_state',
    'particle_input_gc_to_state', 'particle_input_ml_to_state',
    'particle_input_p_to_state', 'particle_input_to_state',
    'particle_ml', 'particle_ml_to_state', 'particle_queue',
    'particle_simd_fo', 'particle_simd_gc', 'particle_simd_ml',
    'particle_state', 'particle_state_to_fo', 'particle_state_to_gc',
    'particle_state_to_ml', 'particle_to_fo_dummy',
    'particle_to_gc_dummy', 'particle_to_ml_dummy', 'plasma_1DS_data',
    'plasma_1DS_eval_dens', 'plasma_1DS_eval_densandtemp',
    'plasma_1DS_eval_temp', 'plasma_1DS_free_offload',
    'plasma_1DS_init', 'plasma_1DS_init_offload',
    'plasma_1DS_offload_data', 'plasma_1D_data',
    'plasma_1D_eval_dens', 'plasma_1D_eval_densandtemp',
    'plasma_1D_eval_temp', 'plasma_1D_free_offload', 'plasma_1D_init',
    'plasma_1D_init_offload', 'plasma_1D_offload_data',
    'plasma_1Dt_data', 'plasma_1Dt_eval_dens',
    'plasma_1Dt_eval_densandtemp', 'plasma_1Dt_eval_temp',
    'plasma_1Dt_free_offload', 'plasma_1Dt_init',
    'plasma_1Dt_init_offload', 'plasma_1Dt_offload_data',
    'plasma_data', 'plasma_eval_dens', 'plasma_eval_densandtemp',
    'plasma_eval_temp', 'plasma_free_offload', 'plasma_get_n_species',
    'plasma_get_species_anum', 'plasma_get_species_charge',
    'plasma_get_species_mass', 'plasma_get_species_znum',
    'plasma_init', 'plasma_init_offload', 'plasma_offload_data',
    'plasma_type', 'plasma_type_1D', 'plasma_type_1DS',
    'plasma_type_1Dt', 'prepare_markers', 'print_marker_summary',
    'real', 'sigma_CX', 'sigma_ioniz', 'sigma_recomb', 'sigmav_BMS',
    'sigmav_CX', 'sigmav_ioniz', 'sigmav_recomb', 'sigmaveff_CX',
    'sigmaveff_ioniz', 'sigmaveff_recomb', 'sim_data', 'sim_init',
    'sim_offload_data', 'simulate', 'simulate_init_offload',
    'simulate_mode_fo', 'simulate_mode_gc', 'simulate_mode_hybrid',
    'simulate_mode_ml', 'size_t', 'struct_c__SA_B_2DS_data',
    'struct_c__SA_B_2DS_offload_data', 'struct_c__SA_B_3DS_data',
    'struct_c__SA_B_3DS_offload_data', 'struct_c__SA_B_GS_data',
    'struct_c__SA_B_GS_offload_data', 'struct_c__SA_B_STS_data',
    'struct_c__SA_B_STS_offload_data', 'struct_c__SA_B_TC_data',
    'struct_c__SA_B_TC_offload_data', 'struct_c__SA_B_field_data',
    'struct_c__SA_B_field_offload_data', 'struct_c__SA_E_1DS_data',
    'struct_c__SA_E_1DS_offload_data', 'struct_c__SA_E_TC_data',
    'struct_c__SA_E_TC_offload_data', 'struct_c__SA_E_field_data',
    'struct_c__SA_E_field_offload_data', 'struct_c__SA_N0_1D_data',
    'struct_c__SA_N0_1D_offload_data', 'struct_c__SA_N0_3D_data',
    'struct_c__SA_N0_3D_offload_data', 'struct_c__SA_afsi_data',
    'struct_c__SA_afsi_thermal_data', 'struct_c__SA_asigma_data',
    'struct_c__SA_asigma_loc_data',
    'struct_c__SA_asigma_loc_offload_data',
    'struct_c__SA_asigma_offload_data', 'struct_c__SA_boozer_data',
    'struct_c__SA_boozer_offload_data', 'struct_c__SA_diag_data',
    'struct_c__SA_diag_offload_data', 'struct_c__SA_diag_orb_data',
    'struct_c__SA_diag_orb_offload_data',
    'struct_c__SA_diag_transcoef_data',
    'struct_c__SA_diag_transcoef_offload_data',
    'struct_c__SA_dist_5D_data', 'struct_c__SA_dist_5D_offload_data',
    'struct_c__SA_dist_6D_data', 'struct_c__SA_dist_6D_offload_data',
    'struct_c__SA_dist_COM_data',
    'struct_c__SA_dist_COM_offload_data',
    'struct_c__SA_dist_rho5D_data',
    'struct_c__SA_dist_rho5D_offload_data',
    'struct_c__SA_dist_rho6D_data',
    'struct_c__SA_dist_rho6D_offload_data',
    'struct_c__SA_input_particle', 'struct_c__SA_interp1D_data',
    'struct_c__SA_interp2D_data', 'struct_c__SA_interp3D_data',
    'struct_c__SA_linint1D_data', 'struct_c__SA_linint3D_data',
    'struct_c__SA_mccc_data', 'struct_c__SA_mhd_data',
    'struct_c__SA_mhd_nonstat_data',
    'struct_c__SA_mhd_nonstat_offload_data',
    'struct_c__SA_mhd_offload_data', 'struct_c__SA_mhd_stat_data',
    'struct_c__SA_mhd_stat_offload_data', 'struct_c__SA_nbi_data',
    'struct_c__SA_nbi_injector', 'struct_c__SA_nbi_offload_data',
    'struct_c__SA_neutral_data', 'struct_c__SA_neutral_offload_data',
    'struct_c__SA_offload_package', 'struct_c__SA_particle',
    'struct_c__SA_particle_gc', 'struct_c__SA_particle_ml',
    'struct_c__SA_particle_queue', 'struct_c__SA_particle_simd_fo',
    'struct_c__SA_particle_simd_gc', 'struct_c__SA_particle_simd_ml',
    'struct_c__SA_particle_state', 'struct_c__SA_plasma_1DS_data',
    'struct_c__SA_plasma_1DS_offload_data',
    'struct_c__SA_plasma_1D_data',
    'struct_c__SA_plasma_1D_offload_data',
    'struct_c__SA_plasma_1Dt_data',
    'struct_c__SA_plasma_1Dt_offload_data',
    'struct_c__SA_plasma_data', 'struct_c__SA_plasma_offload_data',
    'struct_c__SA_sim_data', 'struct_c__SA_sim_offload_data',
    'struct_c__SA_wall_2d_data', 'struct_c__SA_wall_2d_offload_data',
    'struct_c__SA_wall_3d_data', 'struct_c__SA_wall_3d_offload_data',
    'struct_c__SA_wall_data', 'struct_c__SA_wall_offload_data',
    'struct_diag_transcoef_link', 'union_c__SA_input_particle_0',
    'wall_2d_data', 'wall_2d_find_intersection',
    'wall_2d_free_offload', 'wall_2d_hit_wall', 'wall_2d_init',
    'wall_2d_init_offload', 'wall_2d_inside', 'wall_2d_offload_data',
    'wall_3d_data', 'wall_3d_free_offload', 'wall_3d_hit_wall',
    'wall_3d_hit_wall_full', 'wall_3d_init', 'wall_3d_init_octree',
    'wall_3d_init_offload', 'wall_3d_init_tree',
    'wall_3d_offload_data', 'wall_3d_quad_collision',
    'wall_3d_tri_collision', 'wall_3d_tri_in_cube', 'wall_data',
    'wall_free_offload', 'wall_get_n_elements', 'wall_hit_wall',
    'wall_init', 'wall_init_offload', 'wall_offload_data',
    'wall_type', 'wall_type_2D', 'wall_type_3D', 'write_output',
    'write_rungroup']
