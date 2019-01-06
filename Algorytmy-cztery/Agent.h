#ifndef AGENT_H_INCLUDED
#define AGENT_H_INCLUDED

class Agent {
public:
	Agent();
	Agent(std::string);

	std::string getState();

	void setState(std::string);
	int compare(Agent* b);
private:
	std::string state;
};

#endif // AGENT_H_INCLUDED
