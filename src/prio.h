/* > prio.h
 * Headers for printer I/O routines for Survex printer drivers
 * Copyright (C) 1993,1996 Olly Betts
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the Free
 * Software Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <stdio.h>

#include "useful.h"

typedef struct {
   size_t len;
   char *str;
} pstr;

extern void prio_open(const char *fnmPrn);
extern void prio_close(void);
extern void prio_putc(int ch);
extern void prio_printf(const char *format, ...);
extern void prio_print(const char *str);
extern void prio_putpstr(const pstr *p);
extern void prio_putbuf(const void *buf, size_t len);
