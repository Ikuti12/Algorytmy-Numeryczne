#include <iostream>
#include "Agent.h"

Agent::Agent() {}
Agent::Agent(std::string state)
{
    this->state = state;
}

std::string Agent::getState()
{
    return state;
}

void Agent::setState(std::string state)
{
    this->state = state;
}

int Agent::compare(Agent* b)
{
    if (this->getState() == b->getState()) {
        return 0;
    }
    if (this->getState() == "Y") {
        if (b->getState() == "U") {
            b->setState("Y");
            return 1;
        }
        else if (b->getState() == "N") {
            this->setState("U");
            b->setState("U");
            //cout << "WAZNE" << b->getState() << this->getState() << endl;
            return 2;
        }
        else if (b->getState() == "Y") {
            return 0;
        }
    }
    else if (this->getState() == "N") {
        if (b->getState() == "U") {
            b->setState("N");
            return 3;
        }
        else if (b->getState() == "Y") {
            this->setState("U");
            b->setState("U");
            //cout << "WAZNE" << b->getState() << this->getState() << endl;
            return 2;
        }
        else if (b->getState() == "N") {
            return 0;
        }
    }
    else if (b->getState() == "Y") {
        if (this->getState() == "U") {
            this->setState("Y");
            return 1;
        }
        else if (this->getState() == "N") {
            b->setState("U");
            this->setState("U");
            return 2;
        }
        else if (this->getState() == "Y") {
            return 0;
        }
    }

    else if (b->getState() == "N") {
        if (this->getState() == "U") {
            this->setState("N");
            return 3;
        }
        else if (this->getState() == "Y") {
            this->setState("U");
            b->setState("U");
            return 2;
        }
        else if (this->getState() == "N") {
            return 0;
        }
    }
    return 0;
}
