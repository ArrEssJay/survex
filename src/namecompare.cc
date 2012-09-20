/* namecompare.cc */
/* Ordering function for station names */
/* Copyright (C) 1991-2002,2004,2012 Olly Betts
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
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include "namecompare.h"

inline bool u_digit(unsigned ch) {
    return ch < 256 && isdigit((unsigned char)ch);
}

int name_cmp(const wxString &a, const wxString &b, int separator) {
   size_t i = 0;
   while (1) {
      int cha = (i == a.size() ? wxUniChar(0) : a[i]);
      int chb = (i == b.size() ? wxUniChar(0) : b[i]);

      /* done if end of either string */
      if (!cha || !chb) return cha - chb;

      /* check for end of non-numeric prefix */
      if (u_digit(cha)) {
	 /* sort numbers numerically and before non-numbers */
	 size_t sa, sb, ea, eb;
	 int res;

	 if (!u_digit(chb)) return chb == separator ? 1 : -1;

	 sa = i;
	 while (sa != a.size() && a[sa] == '0') sa++;
	 ea = sa;
	 while (ea != a.size() && u_digit((unsigned char)a[ea])) ea++;

	 sb = i;
	 while (sb != b.size() && b[sb] == '0') sb++;
	 eb = sb;
	 while (eb != b.size() && u_digit((unsigned char)b[eb])) eb++;

	 /* shorter sorts first */
	 res = (ea - sa) - (eb - sb);
	 /* same length, all digits, so character value compare sorts
	  * numerically */
	 for (size_t j = sa; !res && j != ea; ++j) {
	    res = a[j] - b[j - sa + sb];
	 }
	 /* more leading zeros sorts first */
	 if (!res) res = sb - sa;
	 if (res) return res;

	 /* if numbers match, sort by suffix */
	 i = ea;
	 continue;
      }

      if (cha != chb) {
	 if (cha == separator) return -1;
	 if (u_digit(chb) || chb == separator) return 1;
	 return cha - chb;
      }

      i++;
   }
}
