/* > listpos.c
 * SURVEX Cave surveying software: stuff to do with stn position output
 * Copyright (C) 1991-1995 Olly Betts
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

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#define PRINT_STN_POS_LIST 1
#define NODESTAT 1

#include "cavern.h"
#include "datain.h"
#include "filename.h"
#include "message.h"
#include "filelist.h"
#include "netbits.h"
#include "listpos.h"
#include "out.h"

static void p_pos(prefix *);

/* Traverse prefix tree depth first starting at from, and
 * calling function fn at each node */
static void
traverse_prefix_tree(prefix *from, void (*fn)(prefix *))
{
   prefix *p;

   fn(from);

   p = from->down;
   if (!p) return;
     
   while (1) {
      fn(p);
      if (p->down) {
	 p = p->down;
      } else {
	 if (!p->right) {
	    do {
	       p = p->up;
	       if (p == from) return; /* got back to start */
	    } while (!p->right);
	 }
	 p = p->right;
      }
   }
}

static FILE *fhPosList;

extern void
list_pos(prefix *from)
{
   fhPosList = safe_fopen_with_ext(fnm_output_base, EXT_SVX_POS, "w");

   /* Output headings line */
   fputsnl(msg(/*( Easting, Northing, Altitude )*/195), fhPosList);

   traverse_prefix_tree(from, p_pos);

   fclose(fhPosList);
}

static void
p_pos(prefix *p)
{
#if PRINT_STN_POS_LIST
   /* check it's not just a stem node, and that the station is fixed */
   if (p->pos && pfx_fixed(p)) {
# ifdef PRINT_STN_ORDER
      /* NB this patch from Leandro Dybal Bertoni <LEANDRO@trieste.fapesp.br>
       * will confuse diffpos a lot so is off by default */
      fprintf(fhPosList, "%2d (%8.2f, %8.2f, %8.2f ) ",
	      p->pos->shape, p->pos->p[0], p->pos->p[1], p->pos->p[2]);
# else
      fprintf(fhPosList, "(%8.2f, %8.2f, %8.2f ) ",
	      p->pos->p[0], p->pos->p[1], p->pos->p[2]);
# endif
      fprint_prefix(fhPosList, p);
      fputnl(fhPosList);
   }
#endif
}

#if NODESTAT
static int *cOrder;
static int icOrderMac;

static void
node_stat(prefix *p)
{
   if (p->pos && pfx_fixed(p)) {
      int order;
      order = p->pos->shape;
      if (!order) warning(/*Unused fixed point '%s'*/73, sprint_prefix(p));
      if (order >= icOrderMac) {
	 int c = order * 2;
	 cOrder = osrealloc(cOrder, c * ossizeof(int));
	 while (icOrderMac < c) cOrder[icOrderMac++] = 0;
      }
      cOrder[order]++;
   }
}

static void
node_stats(prefix *from)
{
   icOrderMac = 8; /* Start value */
   cOrder = osmalloc(icOrderMac * ossizeof(int));
   memset(cOrder, 0, icOrderMac * ossizeof(int));
   traverse_prefix_tree(from, node_stat);
}
#endif

extern void
print_node_stats(FILE *fh)
{
#if NODESTAT
   int c;
   node_stats(root);
   for (c = 0; c < icOrderMac; c++) {
      if (cOrder[c]>0) {
	 char buf[256];
	 sprintf(buf, "%4d %d-%s.", cOrder[c], c,
		 msg(cOrder[c] == 1 ? /*node*/176 : /*nodes*/177));
	 out_puts(buf);
	 fputsnl(buf, fh);
      }
   }
#endif
}
