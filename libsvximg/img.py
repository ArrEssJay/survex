# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.11
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_img', [dirname(__file__)])
        except ImportError:
            import _img
            return _img
        if fp is not None:
            try:
                _mod = imp.load_module('_img', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _img = swig_import_helper()
    del swig_import_helper
else:
    import _img
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0



def img_read_item(*args):
  return _img.img_read_item(*args)
img_read_item = _img.img_read_item
IMG_API_VERSION = _img.IMG_API_VERSION
img_BAD = _img.img_BAD
img_STOP = _img.img_STOP
img_MOVE = _img.img_MOVE
img_LINE = _img.img_LINE
img_LABEL = _img.img_LABEL
img_XSECT = _img.img_XSECT
img_XSECT_END = _img.img_XSECT_END
img_ERROR_INFO = _img.img_ERROR_INFO
img_FLAG_SURFACE = _img.img_FLAG_SURFACE
img_FLAG_DUPLICATE = _img.img_FLAG_DUPLICATE
img_FLAG_SPLAY = _img.img_FLAG_SPLAY
img_SFLAG_SURFACE = _img.img_SFLAG_SURFACE
img_SFLAG_UNDERGROUND = _img.img_SFLAG_UNDERGROUND
img_SFLAG_ENTRANCE = _img.img_SFLAG_ENTRANCE
img_SFLAG_EXPORTED = _img.img_SFLAG_EXPORTED
img_SFLAG_FIXED = _img.img_SFLAG_FIXED
img_SFLAG_ANON = _img.img_SFLAG_ANON
img_SFLAG_WALL = _img.img_SFLAG_WALL
img_FFLAG_EXTENDED = _img.img_FFLAG_EXTENDED
img_XFLAG_END = _img.img_XFLAG_END
img_STYLE_UNKNOWN = _img.img_STYLE_UNKNOWN
img_STYLE_NORMAL = _img.img_STYLE_NORMAL
img_STYLE_DIVING = _img.img_STYLE_DIVING
img_STYLE_CARTESIAN = _img.img_STYLE_CARTESIAN
img_STYLE_CYLPOLAR = _img.img_STYLE_CYLPOLAR
img_STYLE_NOSURVEY = _img.img_STYLE_NOSURVEY
class img_point(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, img_point, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, img_point, name)
    __repr__ = _swig_repr
    __swig_setmethods__["x"] = _img.img_point_x_set
    __swig_getmethods__["x"] = _img.img_point_x_get
    if _newclass:x = _swig_property(_img.img_point_x_get, _img.img_point_x_set)
    __swig_setmethods__["y"] = _img.img_point_y_set
    __swig_getmethods__["y"] = _img.img_point_y_get
    if _newclass:y = _swig_property(_img.img_point_y_get, _img.img_point_y_set)
    __swig_setmethods__["z"] = _img.img_point_z_set
    __swig_getmethods__["z"] = _img.img_point_z_get
    if _newclass:z = _swig_property(_img.img_point_z_get, _img.img_point_z_set)
    def __init__(self): 
        this = _img.new_img_point()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _img.delete_img_point
    __del__ = lambda self : None;
img_point_swigregister = _img.img_point_swigregister
img_point_swigregister(img_point)

class img(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, img, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, img, name)
    __repr__ = _swig_repr
    __swig_setmethods__["label"] = _img.img_label_set
    __swig_getmethods__["label"] = _img.img_label_get
    if _newclass:label = _swig_property(_img.img_label_get, _img.img_label_set)
    __swig_setmethods__["flags"] = _img.img_flags_set
    __swig_getmethods__["flags"] = _img.img_flags_get
    if _newclass:flags = _swig_property(_img.img_flags_get, _img.img_flags_set)
    __swig_setmethods__["title"] = _img.img_title_set
    __swig_getmethods__["title"] = _img.img_title_get
    if _newclass:title = _swig_property(_img.img_title_get, _img.img_title_set)
    __swig_setmethods__["datestamp"] = _img.img_datestamp_set
    __swig_getmethods__["datestamp"] = _img.img_datestamp_get
    if _newclass:datestamp = _swig_property(_img.img_datestamp_get, _img.img_datestamp_set)
    __swig_setmethods__["separator"] = _img.img_separator_set
    __swig_getmethods__["separator"] = _img.img_separator_get
    if _newclass:separator = _swig_property(_img.img_separator_get, _img.img_separator_set)
    __swig_setmethods__["date1"] = _img.img_date1_set
    __swig_getmethods__["date1"] = _img.img_date1_get
    if _newclass:date1 = _swig_property(_img.img_date1_get, _img.img_date1_set)
    __swig_setmethods__["date2"] = _img.img_date2_set
    __swig_getmethods__["date2"] = _img.img_date2_get
    if _newclass:date2 = _swig_property(_img.img_date2_get, _img.img_date2_set)
    __swig_setmethods__["l"] = _img.img_l_set
    __swig_getmethods__["l"] = _img.img_l_get
    if _newclass:l = _swig_property(_img.img_l_get, _img.img_l_set)
    __swig_setmethods__["r"] = _img.img_r_set
    __swig_getmethods__["r"] = _img.img_r_get
    if _newclass:r = _swig_property(_img.img_r_get, _img.img_r_set)
    __swig_setmethods__["u"] = _img.img_u_set
    __swig_getmethods__["u"] = _img.img_u_get
    if _newclass:u = _swig_property(_img.img_u_get, _img.img_u_set)
    __swig_setmethods__["d"] = _img.img_d_set
    __swig_getmethods__["d"] = _img.img_d_get
    if _newclass:d = _swig_property(_img.img_d_get, _img.img_d_set)
    __swig_setmethods__["n_legs"] = _img.img_n_legs_set
    __swig_getmethods__["n_legs"] = _img.img_n_legs_get
    if _newclass:n_legs = _swig_property(_img.img_n_legs_get, _img.img_n_legs_set)
    __swig_setmethods__["length"] = _img.img_length_set
    __swig_getmethods__["length"] = _img.img_length_get
    if _newclass:length = _swig_property(_img.img_length_get, _img.img_length_set)
    __swig_setmethods__["E"] = _img.img_E_set
    __swig_getmethods__["E"] = _img.img_E_get
    if _newclass:E = _swig_property(_img.img_E_get, _img.img_E_set)
    __swig_setmethods__["H"] = _img.img_H_set
    __swig_getmethods__["H"] = _img.img_H_get
    if _newclass:H = _swig_property(_img.img_H_get, _img.img_H_set)
    __swig_setmethods__["V"] = _img.img_V_set
    __swig_getmethods__["V"] = _img.img_V_get
    if _newclass:V = _swig_property(_img.img_V_get, _img.img_V_set)
    __swig_setmethods__["filename_opened"] = _img.img_filename_opened_set
    __swig_getmethods__["filename_opened"] = _img.img_filename_opened_get
    if _newclass:filename_opened = _swig_property(_img.img_filename_opened_get, _img.img_filename_opened_set)
    __swig_setmethods__["is_extended_elevation"] = _img.img_is_extended_elevation_set
    __swig_getmethods__["is_extended_elevation"] = _img.img_is_extended_elevation_get
    if _newclass:is_extended_elevation = _swig_property(_img.img_is_extended_elevation_get, _img.img_is_extended_elevation_set)
    __swig_setmethods__["style"] = _img.img_style_set
    __swig_getmethods__["style"] = _img.img_style_get
    if _newclass:style = _swig_property(_img.img_style_get, _img.img_style_set)
    __swig_setmethods__["fh"] = _img.img_fh_set
    __swig_getmethods__["fh"] = _img.img_fh_get
    if _newclass:fh = _swig_property(_img.img_fh_get, _img.img_fh_set)
    __swig_setmethods__["label_buf"] = _img.img_label_buf_set
    __swig_getmethods__["label_buf"] = _img.img_label_buf_get
    if _newclass:label_buf = _swig_property(_img.img_label_buf_get, _img.img_label_buf_set)
    __swig_setmethods__["buf_len"] = _img.img_buf_len_set
    __swig_getmethods__["buf_len"] = _img.img_buf_len_get
    if _newclass:buf_len = _swig_property(_img.img_buf_len_get, _img.img_buf_len_set)
    __swig_setmethods__["label_len"] = _img.img_label_len_set
    __swig_getmethods__["label_len"] = _img.img_label_len_get
    if _newclass:label_len = _swig_property(_img.img_label_len_get, _img.img_label_len_set)
    __swig_setmethods__["fRead"] = _img.img_fRead_set
    __swig_getmethods__["fRead"] = _img.img_fRead_get
    if _newclass:fRead = _swig_property(_img.img_fRead_get, _img.img_fRead_set)
    __swig_setmethods__["start"] = _img.img_start_set
    __swig_getmethods__["start"] = _img.img_start_get
    if _newclass:start = _swig_property(_img.img_start_get, _img.img_start_set)
    __swig_setmethods__["version"] = _img.img_version_set
    __swig_getmethods__["version"] = _img.img_version_get
    if _newclass:version = _swig_property(_img.img_version_get, _img.img_version_set)
    __swig_setmethods__["survey"] = _img.img_survey_set
    __swig_getmethods__["survey"] = _img.img_survey_get
    if _newclass:survey = _swig_property(_img.img_survey_get, _img.img_survey_set)
    __swig_setmethods__["survey_len"] = _img.img_survey_len_set
    __swig_getmethods__["survey_len"] = _img.img_survey_len_get
    if _newclass:survey_len = _swig_property(_img.img_survey_len_get, _img.img_survey_len_set)
    __swig_setmethods__["pending"] = _img.img_pending_set
    __swig_getmethods__["pending"] = _img.img_pending_get
    if _newclass:pending = _swig_property(_img.img_pending_get, _img.img_pending_set)
    __swig_setmethods__["mv"] = _img.img_mv_set
    __swig_getmethods__["mv"] = _img.img_mv_get
    if _newclass:mv = _swig_property(_img.img_mv_get, _img.img_mv_set)
    __swig_setmethods__["olddate1"] = _img.img_olddate1_set
    __swig_getmethods__["olddate1"] = _img.img_olddate1_get
    if _newclass:olddate1 = _swig_property(_img.img_olddate1_get, _img.img_olddate1_set)
    __swig_setmethods__["olddate2"] = _img.img_olddate2_set
    __swig_getmethods__["olddate2"] = _img.img_olddate2_get
    if _newclass:olddate2 = _swig_property(_img.img_olddate2_get, _img.img_olddate2_set)
    __swig_setmethods__["oldstyle"] = _img.img_oldstyle_set
    __swig_getmethods__["oldstyle"] = _img.img_oldstyle_get
    if _newclass:oldstyle = _swig_property(_img.img_oldstyle_get, _img.img_oldstyle_set)
    def __init__(self): 
        this = _img.new_img()
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _img.delete_img
    __del__ = lambda self : None;
img_swigregister = _img.img_swigregister
img_swigregister(img)

IMG_VERSION_MIN = _img.IMG_VERSION_MIN
IMG_VERSION_MAX = _img.IMG_VERSION_MAX

def img_open_survey(*args):
  return _img.img_open_survey(*args)
img_open_survey = _img.img_open_survey

def img_open_write(*args):
  return _img.img_open_write(*args)
img_open_write = _img.img_open_write

def img_write_item(*args):
  return _img.img_write_item(*args)
img_write_item = _img.img_write_item

def img_write_errors(*args):
  return _img.img_write_errors(*args)
img_write_errors = _img.img_write_errors

def img_rewind(*args):
  return _img.img_rewind(*args)
img_rewind = _img.img_rewind

def img_close(*args):
  return _img.img_close(*args)
img_close = _img.img_close
IMG_NONE = _img.IMG_NONE
IMG_FILENOTFOUND = _img.IMG_FILENOTFOUND
IMG_OUTOFMEMORY = _img.IMG_OUTOFMEMORY
IMG_CANTOPENOUT = _img.IMG_CANTOPENOUT
IMG_BADFORMAT = _img.IMG_BADFORMAT
IMG_DIRECTORY = _img.IMG_DIRECTORY
IMG_READERROR = _img.IMG_READERROR
IMG_WRITEERROR = _img.IMG_WRITEERROR
IMG_TOONEW = _img.IMG_TOONEW

def img_error():
  return _img.img_error()
img_error = _img.img_error
# This file is compatible with both classic and new-style classes.

cvar = _img.cvar

