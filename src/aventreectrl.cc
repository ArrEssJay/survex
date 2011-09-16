//
//  aventreectrl.cc
//
//  Tree control used for the survey tree.
//
//  Copyright (C) 2001, Mark R. Shinwell.
//  Copyright (C) 2001-2003,2005,2006 Olly Betts
//  Copyright (C) 2005 Martin Green
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
//  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
//

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "aventreectrl.h"
#include "mainfrm.h"

BEGIN_EVENT_TABLE(AvenTreeCtrl, wxTreeCtrl)
    EVT_MOTION(AvenTreeCtrl::OnMouseMove)
    EVT_LEAVE_WINDOW(AvenTreeCtrl::OnLeaveWindow)
    EVT_TREE_SEL_CHANGED(-1, AvenTreeCtrl::OnSelChanged)
    EVT_TREE_ITEM_ACTIVATED(-1, AvenTreeCtrl::OnItemActivated)
    EVT_CHAR(AvenTreeCtrl::OnKeyPress)
END_EVENT_TABLE()

AvenTreeCtrl::AvenTreeCtrl(MainFrm* parent, wxWindow* window_parent) :
    wxTreeCtrl(window_parent, -1),
    m_Parent(parent),
    m_Enabled(false),
    m_LastItem(),
    m_BackgroundColour(),
    m_SelValid(false)
{
}

#define TREE_MASK (wxTREE_HITTEST_ONITEMLABEL | wxTREE_HITTEST_ONITEMRIGHT)

void AvenTreeCtrl::OnMouseMove(wxMouseEvent& event)
{
    if (m_Enabled) {
	int flags;
	wxTreeItemId pos = HitTest(event.GetPosition(), flags);
	if (!(flags & TREE_MASK)) {
	    pos = wxTreeItemId();
	}
	if (pos == m_LastItem) return;
	if (pos.IsOk()) {
	    m_Parent->DisplayTreeInfo(GetItemData(pos));
	} else {
	    m_Parent->DisplayTreeInfo(NULL);
	}
    }
}

void AvenTreeCtrl::SetHere(wxTreeItemId pos)
{
    if (pos == m_LastItem) return;

    if (m_LastItem.IsOk()) {
	SetItemBackgroundColour(m_LastItem, m_BackgroundColour);
    }
    if (pos.IsOk()) {
	m_BackgroundColour = GetItemBackgroundColour(pos);
	SetItemBackgroundColour(pos, wxColour(180, 180, 180));
    }
    m_LastItem = pos;
}

void AvenTreeCtrl::OnLeaveWindow(wxMouseEvent&)
{
    if (m_LastItem.IsOk()) {
	SetItemBackgroundColour(m_LastItem, m_BackgroundColour);
	m_LastItem = wxTreeItemId();
    }
    m_Parent->DisplayTreeInfo(NULL);
}

void AvenTreeCtrl::SetEnabled(bool enabled)
{
    m_Enabled = enabled;
}

void AvenTreeCtrl::OnSelChanged(wxTreeEvent& e)
{
    if (m_Enabled) {
	m_Parent->TreeItemSelected(GetItemData(e.GetItem()), false);
    }

    m_SelValid = true;
}

void AvenTreeCtrl::OnItemActivated(wxTreeEvent& e)
{
    if (m_Enabled) {
	m_Parent->TreeItemSelected(GetItemData(e.GetItem()), true);
	// Need to skip to allow double-clicking to work on wxMSW >= 2.8.11.
	e.Skip();
    }
}

bool AvenTreeCtrl::GetSelectionData(wxTreeItemData** data) const
{
    assert(m_Enabled);
    assert(data);

    if (!m_SelValid) {
	return false;
    }

    wxTreeItemId id = GetSelection();
    if (id.IsOk()) {
	*data = GetItemData(id);
    }

    return id.IsOk() && *data;
}

void AvenTreeCtrl::UnselectAll()
{
    m_SelValid = false;
    wxTreeCtrl::UnselectAll();
}

void AvenTreeCtrl::DeleteAllItems()
{
    m_Enabled = false;
    m_LastItem = wxTreeItemId();
    m_SelValid = false;
    wxTreeCtrl::DeleteAllItems();
}

void AvenTreeCtrl::OnKeyPress(wxKeyEvent &e)
{
    switch (e.GetKeyCode()) {
	case WXK_ESCAPE:
	    m_Parent->ClearTreeSelection();
	    break;
	case WXK_RETURN: {
	    wxTreeItemId id = GetSelection();
	    if (id.IsOk()) {
		if (ItemHasChildren(id)) {
		    // If on a branch, expand/contract it.
		    if (IsExpanded(id)) {
			Collapse(id);
		    } else {
			Expand(id);
		    }
		} else {
		    // FIXME if on a station, show information on that station
		    // or something?
		}
	    }
	    break;
	}
	case WXK_LEFT: case WXK_RIGHT: case WXK_UP: case WXK_DOWN:
	case WXK_HOME: case WXK_END: case WXK_PAGEUP: case WXK_PAGEDOWN:
	// On wx 2.6 and earlier, PRIOR/NEXT seem to actually be
	// PAGEUP/PAGEDOWN (on wxGTK at least).  In wx 2.7 and later
	// they're just compatibility aliases, so either they have the
	// same value or aren't defined - either way the code won't
	// compile.
#if !wxCHECK_VERSION(2,7,0)
	case WXK_PRIOR: case WXK_NEXT:
#endif
	    e.Skip();
	    break;
	default:
	    // Pass key event to MainFrm which will pass to GfxCore which will
	    // pass to GUIControl.
	    m_Parent->OnKeyPress(e);
	    break;
    }
}
