/* > readval.h
 * Routines to read a prefix or number from the current input file
 * Copyright (C) 1991-2001 Olly Betts
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
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

prefix *read_prefix_survey(bool fOmit, bool fAllowRoot);
prefix *read_prefix_stn(bool fOmit, bool fAllowRoot);
prefix *read_prefix_stn_check_implicit(bool fOmit, bool fAllowRoot);
real read_numeric(bool fOmit);
unsigned int read_uint(void);
