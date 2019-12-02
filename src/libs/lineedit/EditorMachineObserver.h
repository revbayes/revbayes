/* 
 * File:   EditorMachineObserver.h
 * Author: johadunf
 *
 * Created on May 19, 2014, 8:27 AM
 */

#ifndef EDITORMACHINEOBSERVER_H
#define	EDITORMACHINEOBSERVER_H

#include "EditorState.h"

enum EditorStateChangeType {
    STATE_ADDED, STATE_RELEASED, STATE_CANCELLED, NOOP
};

class EditorMachineObserver {

    public:
    virtual void eventStateChanged(EditorState *state, EditorStateChangeType type) = 0;

};

#endif	/* EDITORMACHINEOBSERVER_H */

