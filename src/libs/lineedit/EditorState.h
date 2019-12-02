#ifndef EDITORSTATE_H
#define	EDITORSTATE_H

#include <string>
#include <vector>
#include <algorithm>

#include "lineeditUtils.h"

#include <boost/assign/std/vector.hpp>

using namespace boost::assign; // bring 'operator+=()' into scope


typedef std::vector<std::string> StringVector;

enum StateType {
    ST_IDLE, ST_DEF_STRING, ST_DEF_LIST, ST_GENERIC_OPERATOR, ST_ACCESSING_MEMBER, ST_DEF_ARGUMENT
};

class EditorState {
public:

    EditorState() {
    }

    std::string getDescription() {
        return description;
    }

    StateType getType() {
        return type;
    }

    std::string getTrigger() {
        return trigger;
    }

    /**
     * Indicates if the cmd input should release this state
     * @param cmd
     * @param type
     * @return 
     */
    virtual bool tryRelease(std::string cmd, StateType type) = 0;

    /**
     * Indicates if the cmd input triggers this state
     * @param cmd
     * @param type
     * @return 
     */
    virtual bool tryHook(std::string cmd, StateType type) = 0;

    StringVector getCompletions() {
        return completions;
    }

    void setCompletions(StringVector completions) {
        this->completions.clear();
        for (unsigned int i = 0; i < completions.size(); i++) {
            this->completions.push_back(completions[i]);
        }
    }

    void addCompletion(std::string completion) {
        completions.push_back(completion);
    }

    void setSubject(std::string subject) {
        this->subject = subject;
    }

    std::string getSubject() const {
        return subject;
    }

protected:

    StringVector completions;
    std::string hookChars; /** any of these triggers this state **/
    StringVector hookStrings;
    char hookChar; // actual trigger
    std::string releaseChars; /** release for this state **/
    std::string subject; /** for example the function name **/
    std::string description;
    std::string trigger;
    StateType type;
};

class StateIdle : public EditorState {
public:

    StateIdle() : EditorState() {
        //        hookChars = "\"";
        //        releaseChars = "\"";
        description = "inIdle";
        type = ST_IDLE;
    }

    virtual bool tryHook(std::string cmd, StateType type) {
        return false; // cannot be triggered
    }

    virtual bool tryRelease(std::string cmd, StateType type) {
        return false; // cannot be released
    }

    virtual bool tryCancel(const char c, StateType type) {
        return false; // cannot be canceled
    }

};

class StateDefiningString : public EditorState {
public:

    StateDefiningString() : EditorState() {
        hookChars = "\"";
        releaseChars = "\"";
        description = "StateDefiningString";
        type = ST_DEF_STRING;
    }

    virtual bool tryHook(std::string cmd, StateType type) {
        if (type == this->type) { // cannot nest strings
            return false;
        }
        if (LineEditUtils().lastCharContains(cmd, hookChars.c_str())) {
            return true;
        }
        return false;
    }

    virtual bool tryRelease(std::string cmd, StateType type) {
        if (type != this->type) { // cannot release unless current state is same type
            return false;
        }
        if (LineEditUtils().lastCharContains(cmd, releaseChars.c_str())) {
            return true;
        }
        return false;
    }
};

class StateDefiningList : public EditorState {
public:

    StateDefiningList() : EditorState() {
        hookChars = "(";
        releaseChars = ")";
        description = "StateDefiningList";
        type = ST_DEF_LIST;
    }

    virtual bool tryHook(std::string cmd, StateType type) {
        // brackets can be nested, but not triggered in middle of a string
        if (type == ST_DEF_STRING) {
            return false;
        }
        if (LineEditUtils().lastCharContains(cmd, hookChars.c_str())) {
            return true;
        }
        return false;
    }

    virtual bool tryRelease(std::string cmd, StateType type) {
        if (!(type == this->type)) { // cannot release unless current state is same type
            return false;
        }
        if (LineEditUtils().lastCharContains(cmd, releaseChars)) {
            return true;
        }
        return false;
    }
};

class StateDefiningArgument : public EditorState {
public:

    StateDefiningArgument() : EditorState() {
        hookChars = ",";
        releaseChars = ",)";
        description = "StateDefiningArgument";
        type = ST_DEF_ARGUMENT;
    }

    virtual bool tryHook(std::string cmd, StateType type) {
        if (type == this->type // cannot nest this state
                || type == ST_DEF_STRING // don't trigger if in middle of a string
                || type != ST_DEF_LIST// don't trigger unless defining a list
                ) {
            return false;
        }

        if (LineEditUtils().lastCharContains(cmd, hookChars.c_str())) {
            return true;
        }
        return false;
    }

    virtual bool tryRelease(std::string cmd, StateType type) {
        if (!(type == this->type
                //|| type == ST_DEF_LIST
                )) { // cannot release unless current state is same type
            return false;
        }

        if (LineEditUtils().lastCharContains(cmd, releaseChars)) {
            return true;
        }
        return false;
    }
};

class StateGenericOperator : public EditorState {
public:

    StateGenericOperator() : EditorState() {
        // important: assign in descending order from longest to shortest strings
        hookStrings += "<-&", "<-", "<<-", ":=", "++", "--", "+=", "-=", "*=", "/=", "&&", "||", "~", "+", "-", "*", "/", "^", "!", "=";
        //releaseChars = " ";
        description = "StateAssigning";
        type = ST_GENERIC_OPERATOR;
    }

    StateGenericOperator(std::string trigger) {
        this->trigger = trigger;
        StateGenericOperator();
    }

    virtual bool tryHook(std::string cmd, StateType type) {
        if (type == this->type // cannot nest this state
                || type == ST_DEF_STRING // don't trigger if in middle of a string
                || type == ST_DEF_LIST// don't trigger if defining a list
                ) {
            return false;
        }
        LineEditUtils::FindLastInfo f = LineEditUtils().findLastOf(hookStrings, cmd, 0);
        if (f.needle != "") {
            trigger = f.needle;
            return true;
        }

        return false;
    }

    virtual bool tryRelease(std::string cmd, StateType type) {
        if (type != this->type) { // cannot release unless current state is same type
            return false;
        }
        // same sequences that triggers also releases
        LineEditUtils::FindLastInfo f = LineEditUtils().findLastOf(hookStrings, cmd, 0);
        if (f.needle != "") {
            trigger = f.needle;
            return true;
        }
        return false;
    }
};

class StateAccessingMember : public EditorState {
public:

    StateAccessingMember() : EditorState() {
        hookChars = ".";
        releaseChars = " "; // space releases this state
        description = "StateAccessingMember";
        type = ST_ACCESSING_MEMBER;
    }

    virtual bool tryHook(std::string cmd, StateType type) {
        if (type == this->type // cannot nest this state
                || type == ST_DEF_STRING // don't trigger if in middle of a string
                ) {
            return false;
        }

        if (LineEditUtils().lastCharContains(cmd, hookChars.c_str())) {
            return true;
        }
        return false;
    }

    virtual bool tryRelease(std::string cmd, StateType type) {
        if (type != this->type) { // cannot release unless current state is same type
            return false;
        }
        if (LineEditUtils().lastCharContains(cmd, releaseChars)) { // space should release this state
            return true;
        }
        return true;
    }
};


#endif	/* EDITORSTATE_H */

