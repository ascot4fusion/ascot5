# -*- coding: utf-8 -*-
#
# TARGET arch is: ['-I/usr/include/hdf5/serial']
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



c_int128 = ctypes.c_ubyte*16
c_uint128 = c_int128
void = None
if ctypes.sizeof(ctypes.c_longdouble) == 16:
    c_long_double_t = ctypes.c_longdouble
else:
    c_long_double_t = ctypes.c_ubyte*16

_libraries = {}
_libraries['libascot.so'] = ctypes.CDLL('./libascot.so')
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
    ('mileage', ctypes.c_double),
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
struct_c__SA_input_particle._fields_ = [
    ('type', input_particle_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('c__SA_input_particle_0', union_c__SA_input_particle_0),
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
class struct_c__SA_B_field_data(Structure):
    pass

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

class struct_c__SA_B_3DS_data(Structure):
    pass

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

struct_c__SA_B_field_data._pack_ = 1 # source:False
struct_c__SA_B_field_data._fields_ = [
    ('type', B_field_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('BGS', struct_c__SA_B_GS_data),
    ('B2DS', struct_c__SA_B_2DS_data),
    ('B3DS', struct_c__SA_B_3DS_data),
    ('BSTS', struct_c__SA_B_STS_data),
    ('BTC', struct_c__SA_B_TC_data),
]

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
a5err = ctypes.c_uint64
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

# values for enumeration 'c__Ea_hdf5_input_options'
c__Ea_hdf5_input_options__enumvalues = {
    1: 'hdf5_input_options',
    2: 'hdf5_input_bfield',
    4: 'hdf5_input_efield',
    8: 'hdf5_input_plasma',
    16: 'hdf5_input_neutral',
    32: 'hdf5_input_wall',
    64: 'hdf5_input_marker',
    128: 'hdf5_input_boozer',
    256: 'hdf5_input_mhd',
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
c__Ea_hdf5_input_options = ctypes.c_uint32 # enum
class struct_c__SA_sim_offload_data(Structure):
    pass

class struct_c__SA_boozer_offload_data(Structure):
    pass

struct_c__SA_boozer_offload_data._pack_ = 1 # source:False
struct_c__SA_boozer_offload_data._fields_ = [
    ('nr', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('r_min', ctypes.c_double),
    ('r_max', ctypes.c_double),
    ('nz', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
    ('z_min', ctypes.c_double),
    ('z_max', ctypes.c_double),
    ('npsi', ctypes.c_int32),
    ('PADDING_2', ctypes.c_ubyte * 4),
    ('psi_min', ctypes.c_double),
    ('psi_max', ctypes.c_double),
    ('psi0', ctypes.c_double),
    ('psi1', ctypes.c_double),
    ('ntheta', ctypes.c_int32),
    ('nthetag', ctypes.c_int32),
    ('r0', ctypes.c_double),
    ('z0', ctypes.c_double),
    ('nrzs', ctypes.c_int32),
    ('offload_array_length', ctypes.c_int32),
]

class struct_c__SA_mhd_offload_data(Structure):
    pass

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


# values for enumeration 'mhd_type'
mhd_type__enumvalues = {
    0: 'mhd_type_stat',
    1: 'mhd_type_nonstat',
}
mhd_type_stat = 0
mhd_type_nonstat = 1
mhd_type = ctypes.c_uint32 # enum
struct_c__SA_mhd_offload_data._pack_ = 1 # source:False
struct_c__SA_mhd_offload_data._fields_ = [
    ('type', mhd_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('stat', struct_c__SA_mhd_stat_offload_data),
    ('nonstat', struct_c__SA_mhd_nonstat_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

class struct_c__SA_E_field_offload_data(Structure):
    pass

class struct_c__SA_E_TC_offload_data(Structure):
    pass

struct_c__SA_E_TC_offload_data._pack_ = 1 # source:False
struct_c__SA_E_TC_offload_data._fields_ = [
    ('Exyz', ctypes.c_double * 3),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_0', ctypes.c_ubyte * 4),
]


# values for enumeration 'E_field_type'
E_field_type__enumvalues = {
    0: 'E_field_type_TC',
    1: 'E_field_type_1DS',
}
E_field_type_TC = 0
E_field_type_1DS = 1
E_field_type = ctypes.c_uint32 # enum
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

struct_c__SA_E_field_offload_data._pack_ = 1 # source:False
struct_c__SA_E_field_offload_data._fields_ = [
    ('type', E_field_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('ETC', struct_c__SA_E_TC_offload_data),
    ('E1DS', struct_c__SA_E_1DS_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

class struct_c__SA_B_field_offload_data(Structure):
    pass

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

struct_c__SA_B_field_offload_data._pack_ = 1 # source:False
struct_c__SA_B_field_offload_data._fields_ = [
    ('type', B_field_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('BGS', struct_c__SA_B_GS_offload_data),
    ('B2DS', struct_c__SA_B_2DS_offload_data),
    ('B3DS', struct_c__SA_B_3DS_offload_data),
    ('BSTS', struct_c__SA_B_STS_offload_data),
    ('BTC', struct_c__SA_B_TC_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

class struct_c__SA_wall_offload_data(Structure):
    pass


# values for enumeration 'wall_type'
wall_type__enumvalues = {
    0: 'wall_type_2D',
    1: 'wall_type_3D',
}
wall_type_2D = 0
wall_type_3D = 1
wall_type = ctypes.c_uint32 # enum
class struct_c__SA_wall_2d_offload_data(Structure):
    pass

struct_c__SA_wall_2d_offload_data._pack_ = 1 # source:False
struct_c__SA_wall_2d_offload_data._fields_ = [
    ('n', ctypes.c_int32),
    ('offload_array_length', ctypes.c_int32),
]

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
    ('PADDING_1', ctypes.c_ubyte * 4),
]

struct_c__SA_wall_offload_data._pack_ = 1 # source:False
struct_c__SA_wall_offload_data._fields_ = [
    ('type', wall_type),
    ('w2d', struct_c__SA_wall_2d_offload_data),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('w3d', struct_c__SA_wall_3d_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

class struct_c__SA_plasma_offload_data(Structure):
    pass


# values for enumeration 'plasma_type'
plasma_type__enumvalues = {
    0: 'plasma_type_1D',
    1: 'plasma_type_1DS',
}
plasma_type_1D = 0
plasma_type_1DS = 1
plasma_type = ctypes.c_uint32 # enum
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

struct_c__SA_plasma_offload_data._pack_ = 1 # source:False
struct_c__SA_plasma_offload_data._fields_ = [
    ('type', plasma_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('plasma_1D', struct_c__SA_plasma_1D_offload_data),
    ('plasma_1DS', struct_c__SA_plasma_1DS_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

class struct_c__SA_diag_offload_data(Structure):
    pass

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

class struct_c__SA_diag_transcoef_offload_data(Structure):
    pass

struct_c__SA_diag_transcoef_offload_data._pack_ = 1 # source:False
struct_c__SA_diag_transcoef_offload_data._fields_ = [
    ('Nmrk', ctypes.c_int64),
    ('Navg', ctypes.c_int32),
    ('recordrho', ctypes.c_int32),
    ('interval', ctypes.c_double),
]

struct_c__SA_diag_offload_data._pack_ = 1 # source:False
struct_c__SA_diag_offload_data._fields_ = [
    ('diagorb_collect', ctypes.c_int32),
    ('dist5D_collect', ctypes.c_int32),
    ('dist6D_collect', ctypes.c_int32),
    ('distrho5D_collect', ctypes.c_int32),
    ('distrho6D_collect', ctypes.c_int32),
    ('diagtrcof_collect', ctypes.c_int32),
    ('diagorb', struct_c__SA_diag_orb_offload_data),
    ('dist5D', struct_c__SA_dist_5D_offload_data),
    ('dist6D', struct_c__SA_dist_6D_offload_data),
    ('distrho5D', struct_c__SA_dist_rho5D_offload_data),
    ('distrho6D', struct_c__SA_dist_rho6D_offload_data),
    ('diagtrcof', struct_c__SA_diag_transcoef_offload_data),
    ('offload_dist5D_index', ctypes.c_int32),
    ('offload_dist6D_index', ctypes.c_int32),
    ('offload_distrho5D_index', ctypes.c_int32),
    ('offload_distrho6D_index', ctypes.c_int32),
    ('offload_diagorb_index', ctypes.c_int32),
    ('offload_diagtrcof_index', ctypes.c_int32),
    ('offload_dist_length', ctypes.c_int32),
    ('offload_array_length', ctypes.c_int32),
]

class struct_c__SA_neutral_offload_data(Structure):
    pass


# values for enumeration 'neutral_type'
neutral_type__enumvalues = {
    0: 'neutral_type_3D',
}
neutral_type_3D = 0
neutral_type = ctypes.c_uint32 # enum
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

struct_c__SA_neutral_offload_data._pack_ = 1 # source:False
struct_c__SA_neutral_offload_data._fields_ = [
    ('type', neutral_type),
    ('PADDING_0', ctypes.c_ubyte * 4),
    ('N03D', struct_c__SA_N0_3D_offload_data),
    ('offload_array_length', ctypes.c_int32),
    ('PADDING_1', ctypes.c_ubyte * 4),
]

struct_c__SA_sim_offload_data._pack_ = 1 # source:False
struct_c__SA_sim_offload_data._fields_ = [
    ('B_offload_data', struct_c__SA_B_field_offload_data),
    ('E_offload_data', struct_c__SA_E_field_offload_data),
    ('plasma_offload_data', struct_c__SA_plasma_offload_data),
    ('neutral_offload_data', struct_c__SA_neutral_offload_data),
    ('wall_offload_data', struct_c__SA_wall_offload_data),
    ('boozer_offload_data', struct_c__SA_boozer_offload_data),
    ('mhd_offload_data', struct_c__SA_mhd_offload_data),
    ('diag_offload_data', struct_c__SA_diag_offload_data),
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
    ('disable_gctransform', ctypes.c_int32),
    ('disable_energyccoll', ctypes.c_int32),
    ('disable_pitchccoll', ctypes.c_int32),
    ('disable_gcdiffccoll', ctypes.c_int32),
    ('endcond_active', ctypes.c_int32),
    ('endcond_max_simtime', ctypes.c_double),
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
    ('PADDING_1', ctypes.c_ubyte * 4),
]

hdf5_interface_read_input = _libraries['libascot.so'].hdf5_interface_read_input
hdf5_interface_read_input.restype = ctypes.c_int32
hdf5_interface_read_input.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.c_int32, ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(struct_c__SA_input_particle)), ctypes.POINTER(ctypes.c_int32)]
hdf5_interface_init_results = _libraries['libascot.so'].hdf5_interface_init_results
hdf5_interface_init_results.restype = ctypes.c_int32
hdf5_interface_init_results.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.POINTER(ctypes.c_char)]
hdf5_interface_write_state = _libraries['libascot.so'].hdf5_interface_write_state
hdf5_interface_write_state.restype = ctypes.c_int32
hdf5_interface_write_state.argtypes = [ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_char), integer, ctypes.POINTER(struct_c__SA_particle_state)]
hdf5_interface_write_diagnostics = _libraries['libascot.so'].hdf5_interface_write_diagnostics
hdf5_interface_write_diagnostics.restype = ctypes.c_int32
hdf5_interface_write_diagnostics.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_char)]
hid_t = ctypes.c_int64
hdf5_get_active_qid = _libraries['libascot.so'].hdf5_get_active_qid
hdf5_get_active_qid.restype = ctypes.c_int32
hdf5_get_active_qid.argtypes = [hid_t, ctypes.POINTER(ctypes.c_char), ctypes.c_char * 11]
hdf5_generate_qid = _libraries['libascot.so'].hdf5_generate_qid
hdf5_generate_qid.restype = None
hdf5_generate_qid.argtypes = [ctypes.POINTER(ctypes.c_char)]
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
class struct_c__SA_offload_package(Structure):
    pass

struct_c__SA_offload_package._pack_ = 1 # source:False
struct_c__SA_offload_package._fields_ = [
    ('offload_array_length', ctypes.c_uint64),
    ('unpack_pos', ctypes.c_uint64),
]

offload = _libraries['libascot.so'].offload
offload.restype = ctypes.c_int32
offload.argtypes = [ctypes.POINTER(struct_c__SA_sim_offload_data), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.c_int32, ctypes.POINTER(ctypes.c_char), ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.POINTER(struct_c__SA_input_particle)), ctypes.POINTER(ctypes.c_int32), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(struct_c__SA_offload_package), ctypes.POINTER(ctypes.POINTER(struct_c__SA_particle_state)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double)), ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
__all__ = \
    ['B_field_type', 'B_field_type_2DS', 'B_field_type_3DS',
    'B_field_type_GS', 'B_field_type_STS', 'B_field_type_TC',
    'E_field_type', 'E_field_type_1DS', 'E_field_type_TC', 'a5err',
    'c__Ea_hdf5_input_options', 'hdf5_generate_qid',
    'hdf5_get_active_qid', 'hdf5_input_bfield', 'hdf5_input_boozer',
    'hdf5_input_efield', 'hdf5_input_marker', 'hdf5_input_mhd',
    'hdf5_input_neutral', 'hdf5_input_options', 'hdf5_input_plasma',
    'hdf5_input_wall', 'hdf5_interface_init_results',
    'hdf5_interface_read_input', 'hdf5_interface_write_diagnostics',
    'hdf5_interface_write_state', 'hid_t', 'input_particle',
    'input_particle_type', 'input_particle_type_gc',
    'input_particle_type_ml', 'input_particle_type_p',
    'input_particle_type_s', 'integer', 'mhd_type',
    'mhd_type_nonstat', 'mhd_type_stat', 'mpi_gather_diag',
    'mpi_gather_particlestate', 'mpi_interface_finalize',
    'mpi_interface_init', 'mpi_my_particles', 'neutral_type',
    'neutral_type_3D', 'offload', 'particle', 'particle_copy_fo',
    'particle_copy_gc', 'particle_copy_ml', 'particle_cycle_fo',
    'particle_cycle_gc', 'particle_cycle_ml', 'particle_fo_to_gc',
    'particle_fo_to_state', 'particle_gc', 'particle_gc_to_state',
    'particle_input_to_state', 'particle_ml', 'particle_ml_to_state',
    'particle_queue', 'particle_simd_fo', 'particle_simd_gc',
    'particle_simd_ml', 'particle_state', 'particle_state_to_fo',
    'particle_state_to_gc', 'particle_state_to_ml',
    'particle_to_fo_dummy', 'particle_to_gc_dummy',
    'particle_to_ml_dummy', 'plasma_type', 'plasma_type_1D',
    'plasma_type_1DS', 'real', 'struct_c__SA_B_2DS_data',
    'struct_c__SA_B_2DS_offload_data', 'struct_c__SA_B_3DS_data',
    'struct_c__SA_B_3DS_offload_data', 'struct_c__SA_B_GS_data',
    'struct_c__SA_B_GS_offload_data', 'struct_c__SA_B_STS_data',
    'struct_c__SA_B_STS_offload_data', 'struct_c__SA_B_TC_data',
    'struct_c__SA_B_TC_offload_data', 'struct_c__SA_B_field_data',
    'struct_c__SA_B_field_offload_data',
    'struct_c__SA_E_1DS_offload_data',
    'struct_c__SA_E_TC_offload_data',
    'struct_c__SA_E_field_offload_data',
    'struct_c__SA_N0_3D_offload_data',
    'struct_c__SA_boozer_offload_data',
    'struct_c__SA_diag_offload_data',
    'struct_c__SA_diag_orb_offload_data',
    'struct_c__SA_diag_transcoef_offload_data',
    'struct_c__SA_dist_5D_offload_data',
    'struct_c__SA_dist_6D_offload_data',
    'struct_c__SA_dist_rho5D_offload_data',
    'struct_c__SA_dist_rho6D_offload_data',
    'struct_c__SA_input_particle', 'struct_c__SA_interp2D_data',
    'struct_c__SA_interp3D_data', 'struct_c__SA_linint1D_data',
    'struct_c__SA_mhd_nonstat_offload_data',
    'struct_c__SA_mhd_offload_data',
    'struct_c__SA_mhd_stat_offload_data',
    'struct_c__SA_neutral_offload_data',
    'struct_c__SA_offload_package', 'struct_c__SA_particle',
    'struct_c__SA_particle_gc', 'struct_c__SA_particle_ml',
    'struct_c__SA_particle_queue', 'struct_c__SA_particle_simd_fo',
    'struct_c__SA_particle_simd_gc', 'struct_c__SA_particle_simd_ml',
    'struct_c__SA_particle_state',
    'struct_c__SA_plasma_1DS_offload_data',
    'struct_c__SA_plasma_1D_offload_data',
    'struct_c__SA_plasma_offload_data',
    'struct_c__SA_sim_offload_data',
    'struct_c__SA_wall_2d_offload_data',
    'struct_c__SA_wall_3d_offload_data',
    'struct_c__SA_wall_offload_data', 'union_c__SA_input_particle_0',
    'wall_type', 'wall_type_2D', 'wall_type_3D']
