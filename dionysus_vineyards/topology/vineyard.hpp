/* Implementations */

#include <fstream>
#include <sstream>
#include <math.h>

#include "utilities/log.h"

#ifdef LOGGING
static rlog::RLogChannel* rlVineyard =          DEF_CHANNEL("topology/vineyard", rlog::Log_Debug);
#endif // LOGGING

template<class I, class It, class E>
void
Vineyard<I,It,E>::
start_vines(Iterator bg, Iterator end)
{
    AssertMsg(evaluator != 0, "Cannot start vines with a null evaluator");
    for (Iterator cur = bg; cur != end; ++cur)
    {
        if (!cur->sign()) continue;
        Dimension dim = evaluator->dimension(cur);
        
        if (dim >= vines.size())
        {
            AssertMsg(dim == vines.size(), "New dimension has to be contiguous");
            vines.push_back(VineList());
            vines_vector.push_back(boost::prior(vines.end()));
        }

        start_vine(cur);
        record_knee(cur);
    }
}

template<class I, class It, class E>
void
Vineyard<I,It,E>::
switched(Index i, Index j)
{
    rLog(rlVineyard, "Switching vines");

    Vine* i_vine = i->vine();
    Vine* j_vine = j->vine();
    i->set_vine(j_vine);
    j->set_vine(i_vine);

    // rLog(rlVineyard, "  %x %x %x %x", i->pair, i->pair->pair, j->pair, j->pair->pair);
    // rLog(rlVineyard, "  %x %x %x %x", i_vine, i->pair->vine(), j_vine, j->pair->vine());

    // Since the pairing has already been updated, the following assertions should be true
    AssertMsg(i->vine() == i->pair->vine(), "i's vine must match the vine of its pair");
    AssertMsg(j->vine() == j->pair->vine(), "j's vine must match the vine of its pair");

    if (!i->sign()) i = i->pair;
    if (!j->sign()) j = j->pair;

    // std::cout << "i sign: " << i->sign() << std::endl;
    // std::cout << "j sign: " << j->sign() << std::endl;

    record_knee(i);
    record_knee(j);
}

template<class I, class It, class E>
template<class Iter>
void
Vineyard<I,It,E>::
start_vine(Iter i)
{
    rLog(rlVineyard, "Starting new vine");
    AssertMsg(i->sign(), "Can only start vines for positive simplices");
        
    Dimension dim = evaluator->dimension(i);
    vines_vector[dim]->push_back(Vine());
    i->set_vine(&vines_vector[dim]->back());
    i->pair->set_vine(i->vine());
}
    
/// Records the current diagram in the vineyard
template<class I, class It, class E>
void 
Vineyard<I,It,E>::
record_diagram(Iterator bg, Iterator end)
{
    rLog(rlVineyard, "Entered: record_diagram()");
    AssertMsg(evaluator != 0, "Cannot record diagram with a null evaluator");
    
    for (Iterator i = bg; i != end; ++i)
    {
        AssertMsg(i->vine() != 0, "Cannot process a null vine in record_diagram");
        if (!i->sign())     continue;
        record_knee(i);
    }
}


template<class I, class It, class E>
void            
Vineyard<I,It,E>::
save_edges(const std::string& filename, bool skip_infinite) const
{
    for (unsigned int i = 0; i < vines_vector.size(); ++i)
    {
        std::ostringstream os; os << i;
        std::string fn = filename + os.str() + ".edg";
        std::ofstream out(fn.c_str());
        for (typename VineList::const_iterator vi = vines_vector[i]->begin(); vi != vines_vector[i]->end(); ++vi)
            for (typename Vine::const_iterator ki = vi->begin(), kiprev = ki++; ki != vi->end(); kiprev = ki++)
            {
                if (skip_infinite && (kiprev->is_infinite() || ki->is_infinite()))
                {
                    std::cerr << "Warning: skipping an infinite knee in save_edges() in dimension " << i << std::endl;
                    continue;
                }
                out << kiprev->birth << ' ' << kiprev->death << ' ' << kiprev->time << std::endl;
                out << ki->birth << ' ' << ki->death << ' ' << ki->time << std::endl;
            }
        out.close();
    }
}

template<class I, class It, class E>
void            
Vineyard<I,It,E>::
save_vines(const std::string& filename, bool skip_infinite) const
{
    for (unsigned int i = 0; i < vines_vector.size(); ++i)
    {
        std::ostringstream os; os << i;
        std::string fn = filename + os.str() + ".vin";
        std::ofstream out(fn.c_str());
        for (typename VineList::const_iterator vi = vines_vector[i]->begin(); vi != vines_vector[i]->end(); ++vi)
        {
            for (typename Vine::const_iterator ki = vi->begin(); ki != vi->end(); ki++)
            {
                if (skip_infinite && ki->is_infinite())
                {
                    std::cerr << "Warning: skipping an infinite knee in save_edges() in dimension " << i << std::endl;
                    continue;
                }
                out << ki->birth << ' ' << ki->death << ' ' << ki->time << " ";
            }
            out << std::endl;
        }
        out.close();
    }
}

template<class I, class It, class E>
std::vector<std::vector<std::vector<double>>>            
Vineyard<I,It,E>::
get_vines(const int& discard_infinite) const
{
    //std::cout << discard_infinite << std::endl;
    std::vector<std::vector<std::vector<double>>> VVV;
    for (unsigned int i = 0; i < vines_vector.size(); ++i)
    {
        std::vector<std::vector<double>> VV;
        for (typename VineList::const_iterator vi = vines_vector[i]->begin(); vi != vines_vector[i]->end(); ++vi)
        {
            std::vector<double> V;
            for (typename Vine::const_iterator ki = vi->begin(); ki != vi->end(); ki++)
            {
                if (discard_infinite && ki->is_infinite())
                {
                    //std::cout << "discarding stuff" << std::endl;
                    continue;
                }
                //std::cout << ki->birth << " " << ki->death << " " << ki->time << std::endl;
                V.push_back(ki->birth); V.push_back(ki->death); V.push_back(ki->time);
            }
            VV.push_back(V);
        }
        VVV.push_back(VV);
    }
    return VVV;
}

template<class I, class It, class E>
std::vector<std::vector<std::vector<double>>>            
Vineyard<I,It,E>::
get_dgms(const int& discard_infinite, const int& num) const
{
    //std::cout << discard_infinite << std::endl;
    std::vector<std::vector<std::vector<double>>> VVV;
    for (unsigned int i = 0; i < vines_vector.size(); ++i)
    {
        std::vector<std::vector<double>> VV(num);
        for (typename VineList::const_iterator vi = vines_vector[i]->begin(); vi != vines_vector[i]->end(); ++vi)
        {
            for (typename Vine::const_iterator ki = vi->begin(); ki != vi->end(); ki++)
            {
                if (discard_infinite && ki->is_infinite())
                {
                    //std::cout << "discarding stuff" << std::endl;
                    continue;
                }
                //std::cout << ki->birth << " " << ki->death << " " << ki->time << std::endl;
                int absulate = abs(ki->time);
                if (absulate == ki->time){VV[ki->time].push_back(ki->birth); VV[ki->time].push_back(ki->death);}
            }
        }
        VVV.push_back(VV);
    }
    return VVV;
}

/// Records a knee for the given simplex
template<class I, class It, class E>
template<class Iter>
void
Vineyard<I,It,E>::
record_knee(Iter i)
{
    rLog(rlVineyard, "Entered record_knee()");
    AssertMsg(evaluator != 0, "Cannot record knee with a null evaluator");
    AssertMsg(i->vine() != 0, "Cannot add a knee to a null vine");
    AssertMsg(i->sign(), "record_knee() must be called on a positive simplex");
    
    if (i->unpaired())
        i->vine()->add((*evaluator)(i), Infinity, evaluator->time());
    else
    {
        rLog(rlVineyard, "Creating knee");
        Knee k((*evaluator)(i), (*evaluator)((i->pair)), evaluator->time());
        rLog(rlVineyard, "Knee created: %s", tostring(k).c_str());
        rLog(rlVineyard, "Vine: %s", tostring(*(i->vine())).c_str());

        if (!k.is_diagonal() || i->vine()->empty())         // non-diagonal k, or empty vine
        {
            rLog(rlVineyard, "Extending a vine");
            i->vine()->add(k);
        }
        else if (i->vine()->back().is_diagonal())           // last knee is diagonal
        {
            AssertMsg(i->vine()->size() == 1, "Only first knee may be diagonal for a live vine");
            rLog(rlVineyard, "Overwriting first diagonal knee");
            i->vine()->back() = k;
        } else                                              // finish this vine
        {
            rLog(rlVineyard, "Finishing a vine");
            i->vine()->add(k);
            start_vine(i);
            i->vine()->add(k);
        }
    }
    
    rLog(rlVineyard, "Leaving record_knee()");
}
