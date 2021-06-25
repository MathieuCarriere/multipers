#include <ctime>
#include <cstdio>

Counter::
Counter(const std::string& full_name,
		CounterType freq):
		count(0), frequency(freq), trigger(this), full_name_(full_name)
{ 
	if (isatty(STDOUT_FILENO)) 
	{
		start_color = green_color; 
		finish_color = normal_color;
	} 
	else
	{	start_color = finish_color = empty_string;	}
}


Counter*
Counter::
get_child(const std::string& path, std::string::size_type pos)
{
	if (pos >= path.size())
		return this;

	std::string::size_type slash_pos = path.find('/', pos);
	if (slash_pos == std::string::npos)
		slash_pos = path.size();

	std::string child_name = path.substr(pos, slash_pos - pos);
	SubCounterMap::iterator child = subcounters_.find(child_name);

	if (child != subcounters_.end())
		return child->second->get_child(path, slash_pos + 1);
	else	
		return (subcounters_[child_name] = new Counter(path.substr(0, slash_pos)))->get_child(path, slash_pos + 1);
}

Counter::
~Counter()
{ 
	if (full_name_ == "" && (!subcounters_.empty() || count))
		print(); 

	for (SubCounterMap::iterator cur = subcounters_.begin(); cur != subcounters_.end(); ++cur)
		delete cur->second;
}

inline
void
Counter::
print()
{
	time_t rawtime; time(&rawtime);
	struct tm* timeinfo = localtime(&rawtime);

	printf("%s(%02i:%02i:%02i)%s ",
		    start_color,	
			timeinfo->tm_hour,
			timeinfo->tm_min,
			timeinfo->tm_sec,
			finish_color);
	std::cout << "Counter [" << full_name_ << "]: " << count << std::endl;
	for (SubCountMap::const_iterator cur = subcount.begin(); cur != subcount.end(); ++cur)
		std::cout << "    " << cur->first << ": " << cur->second << std::endl;
	for (SubCounterMap::iterator cur = subcounters_.begin(); cur != subcounters_.end(); ++cur)
		cur->second->print();
}	

