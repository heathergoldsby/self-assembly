/* main.cpp
 * 
 * This file is part of Self-Assembly.
 * 
 * Copyright 2012 David B. Knoester.
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/max.hpp>

#include <ea/artificial_life/artificial_life.h>
#include <ea/artificial_life/hardware.h>
#include <ea/artificial_life/isa.h>
#include <ea/artificial_life/spatial.h>
#include <ea/artificial_life/datafiles/reactions.h>
#include <ea/artificial_life/datafiles/generation_priority.h>
#include <ea/selection/random.h>
#include <ea/selection/proportionate.h>
#include <ea/meta_population.h>
#include <ea/cmdline_interface.h>

using namespace ea;


LIBEA_MD_DECL(FF_FITNESS, "ea.french_flag.fitness", double);


/*! French flag-based population competition.
 
 blue, white, red verical stripes...
 */
template <typename EA>
struct french_flag : periodic_event<METAPOP_COMPETITION_PERIOD,EA> {
    french_flag(EA& ea) : periodic_event<METAPOP_COMPETITION_PERIOD,EA>(ea), _df("french_flag.dat") {
        _df.add_field("update")
        .add_field("mean_fitness")
        .add_field("max_fitness");
    }
    
    virtual ~french_flag() {
    }
    
    virtual void operator()(EA& ea) {        
        using namespace boost::accumulators;
        accumulator_set<double, stats<tag::mean, tag::max> > fit;
        
        // calculate "fitness":
        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {
            double fitness=0.0;
            for(typename EA::individual_type::population_type::iterator j=(*i)->population().begin(); j!=(*i)->population().end(); ++j) {
                typename EA::individual_type::individual_type& org=**j;
                if(exists<LOCATION_COLOR>(*org.location())) {
                    // what color is it?
                    int color = get<LOCATION_COLOR>(*org.location());
                    int stripe = org.location()->x/(get<SPATIAL_X>(ea)/3); // integer division: 0,1,2
                    
                    if(color == stripe) {
                        ++fitness;
                    }
                }
            }
            fit(fitness);
            put<FF_FITNESS>(pow(1.08,fitness),**i);
        }        
        
        _df.write(ea.current_update())
        .write(mean(fit))
        .write(max(fit))
        .endl();

        // how many survivors?
        std::size_t n = static_cast<std::size_t>((1.0 - get<REPLACEMENT_RATE_P>(ea)) * get<META_POPULATION_SIZE>(ea));
        
        // select individuals for survival:
        typename EA::population_type survivors;
        select_n<selection::random>(ea.population(), survivors, n, ea);
        
        // how many offspring?
        n = get<META_POPULATION_SIZE>(ea) - survivors.size();
        
        // select the parents:
        typename EA::population_type parents;
        select_n<selection::proportionate<attributes::meta_data<FF_FITNESS> > >(survivors, parents, n, ea);
        
        // now, "recombine" the parents to produce offspring:
        typename EA::population_type offspring;
        for(typename EA::population_type::iterator i=parents.begin(); i!=parents.end(); ++i) {
            typename EA::individual_ptr_type p(new typename EA::individual_type());

            // setup the population (really, an ea):
            p->md() = ea.md();
            p->rng().reset(ea.rng()(std::numeric_limits<int>::max()));
            p->initialize();
            
            // grab a copy of the first individual:
            typename EA::individual_type::individual_type germ(**((*i)->population().begin()));
            
            // mutate it:
            mutate(germ,*p);

            // and fill up the offspring population with copies of the germ:
            for(std::size_t j=0; j<get<POPULATION_SIZE>(ea); ++j) {
                typename EA::individual_type::individual_ptr_type o=make_population_entry(germ,*p);
                p->population().push_back(o);
                p->env().insert(o);
            }
            
            offspring.push_back(p);
        }
        
        //reset the survivors!
                
        
        // add the offspring to the list of survivors:
        survivors.insert(survivors.end(), offspring.begin(), offspring.end());
                
        // and swap 'em in for the current population:
        std::swap(ea.population(), survivors);
    }

    datafile _df;
};


template <typename EA>
struct french_flag_movie : public ea::analysis::unary_function<EA> {
    static const char* name() { return "french_flag_movie"; }
    
    virtual void operator()(EA& ea) {
        using namespace ea;
        using namespace ea::analysis;
        using namespace boost::accumulators;

        typename EA::population_type::iterator dominant=ea.population().begin();
        for(typename EA::population_type::iterator i=ea.population().begin(); i!=ea.population().end(); ++i) {
            if(get<FF_FITNESS>(**i) > get<FF_FITNESS>(**dominant)) {
                dominant = i;
            }
        }

        // build a new ea that we transfer the individuals into (because serialization is broken, grr...
        typename EA::individual_ptr_type p(new typename EA::individual_type());
            
        // setup the population (really, an ea):
        p->md() = ea.md();
        p->rng().reset(ea.rng()(std::numeric_limits<int>::max()));
        p->initialize();
            
        // grab a copy of the first individual:        
        typename EA::individual_type::individual_type germ(**((*dominant)->population().begin()));
            
        // fill up the population with copies of the germ:
        for(std::size_t j=0; j<get<POPULATION_SIZE>(ea); ++j) {
            typename EA::individual_type::individual_ptr_type o=make_population_entry(germ,*p);
            p->population().push_back(o);
            p->env().insert(o);
        }
            
        // run the ea for a bit...
        for(unsigned int i=0; i<get<METAPOP_COMPETITION_PERIOD>(ea); ++i) {
            p->update();
        }
    }
};


/*! Artificial life simulation definition.
 */
typedef artificial_life<
hardware, isa, spatial, first_neighbor, round_robin,
mutation::per_site<mutation::uniform_integer>, task_library, organism, population,
alife_population<nopx_ancestor>
> ea_type;

typedef meta_population<
ea_type> mea_type;

/*! 
 */
template <typename EA>
class cli : public cmdline_interface<EA> {
public:
    virtual void configure(EA& ea) {
    }
    
    virtual void gather_options() {
        add_option<SPATIAL_X>(this);
        add_option<SPATIAL_Y>(this);
        add_option<META_POPULATION_SIZE>(this);
        add_option<METAPOP_COMPETITION_PERIOD>(this);
        add_option<REPLACEMENT_RATE_P>(this);
        add_option<POPULATION_SIZE>(this);
        add_option<INITIAL_POPULATION_SIZE>(this);
        add_option<REPRESENTATION_SIZE>(this);
        add_option<SCHEDULER_TIME_SLICE>(this);
        add_option<MUTATION_PER_SITE_P>(this);
        add_option<MUTATION_INSERTION_P>(this);
        add_option<MUTATION_DELETION_P>(this);
        add_option<MUTATION_UNIFORM_INT_MIN>(this);
        add_option<MUTATION_UNIFORM_INT_MAX>(this);
        add_option<RUN_UPDATES>(this);
        add_option<RUN_EPOCHS>(this);
        add_option<CHECKPOINT_PREFIX>(this);        
        add_option<RNG_SEED>(this);
        add_option<RECORDING_PERIOD>(this);
        
        // analysis options
        add_option<ANALYSIS_INPUT>(this);
        add_option<ANALYSIS_OUTPUT>(this);
        add_option<ANALYSIS_ROUNDS>(this);
    }
    
    virtual void gather_tools() {
        add_tool<french_flag_movie>(this);
    }
    
    virtual void gather_events(EA& ea) {
        add_event<french_flag>(this,ea);
//        add_event<datafiles::record_reactions_event>(this,ea);
//        add_event<datafiles::generation_priority>(this,ea);
    };
};
LIBEA_CMDLINE_INSTANCE(mea_type, cli);
