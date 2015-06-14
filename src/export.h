/* export.h
 * Export to CAD-like formats (DXF, Skencil, SVG, EPS, HPGL) and also Compass
 * PLT.
 */

/* Copyright (C) 2004,2005,2012,2014,2015 Olly Betts
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
 */

#ifndef SURVEX_INCLUDED_EXPORT_H
#define SURVEX_INCLUDED_EXPORT_H

#include "wx.h"

#include <time.h>

class MainFrm;

typedef enum {
    FMT_DXF,
    FMT_EPS,
    FMT_GPX,
    FMT_HPGL,
    FMT_JSON,
    FMT_KML,
    FMT_PLT,
    FMT_POS,
    FMT_SK,
    FMT_SVG
} export_format;

#define LEGS		0x00000001
#define SURF		0x00000002
#define STNS		0x00000004
#define LABELS		0x00000008
#define XSECT		0x00000010
#define WALL1		0x00000020
#define WALL2		0x00000040
#define WALLS (WALL1|WALL2)
#define PASG		0x00000080
#define EXPORT_3D	0x10000000
#define CENTRED		0x00000200
#define ENTS		0x00000400
#define FIXES		0x00000800
#define EXPORTS		0x00001000
#define PROJ		0x00002000
#define GRID		0x00004000
#define TEXT_HEIGHT	0x00008000
#define MARKER_SIZE	0x00010000
#define SCALE		0x00020000
#define FULL_COORDS	0x00040000

bool Export(const wxString &fnm_out, const wxString &title,
	    const wxString &datestamp, time_t datestamp_numeric,
	    const MainFrm * mainfrm,
	    double pan, double tilt, int show_mask, export_format format,
	    const char * input_projection,
	    double grid_, double text_height_, double marker_size_,
	    double scale);

#endif
