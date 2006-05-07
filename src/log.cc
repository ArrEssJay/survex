//  log.cc - Error log window for Aven.
//
//  Copyright (C) 2006 Olly Betts
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "aven.h"
#include "log.h"
#include "gla.h"

//#include "message.h"

// FIXME: Make window title available for translation.
MyLogWindow::MyLogWindow()
    : wxLogWindow(NULL, APP_NAME" Error Log", false, false), first(true)
{
}

void MyLogWindow::DoLogString(const wxChar *msg, time_t timestamp) {
    if (first) {
	wxLogWindow::DoLogString(GetGLSystemDescription().c_str(), timestamp);
	first = false;
    }
    wxLogWindow::DoLogString(msg, timestamp);
    Show();
}
