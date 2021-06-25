template<class S>
typename ConeSimplex<S>::Cycle
ConeSimplex<S>::boundary() const
{
	Cycle bdry;
	typename Parent::Cycle pbdry = Parent::boundary();

	for (typename Parent::Cycle::const_iterator cur = pbdry.begin(); cur != pbdry.end(); ++cur)
		bdry.push_back(Self(*cur, coned_));
	
	if (coned_)
		bdry.push_back(Self(*this, false));
	
	return bdry;
}

template<class S>
std::ostream&
ConeSimplex<S>::operator<<(std::ostream& out) const
{
	Parent::operator<<(out) << ' ';
	if (coned_) out << "[coned]";
	return out;
}
