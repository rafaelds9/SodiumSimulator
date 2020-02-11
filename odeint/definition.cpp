/* Boost definition.cpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 The code here is only used for the definitions in the article.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/


template<
    class Container ,
    class Time = double ,
    class Traits = container_traits< Container >
    >
class ode_step
{
    // provide basic typedefs
    //
 public:

    typedef unsigned short order_type;
    typedef Time time_type;
    typedef Traits traits_type;
    typedef typename traits_type::container_type container_type;
    typedef typename traits_type::value_type value_type;

    // public interface
 public:

    // return the order of the stepper
    order_type order_step() const;

    // standard constructor
    ode_step( void );

    // adjust the size of the internal arrays
    void adjust_size( const container_type &x );

    // performs one step
    template< class DynamicalSystem >
    void do_step( DynamicalSystem &system ,
                  container_type &x ,
                  const container_type &dxdt ,
                  time_type t ,
                  time_type dt );

    // performs one step
    template< class DynamicalSystem >
    void do_step( DynamicalSystem &system ,
                  container_type &x ,
                  time_type t ,
                  time_type dt );
};








template<
    class Container ,
    class Time = double ,
    class Traits = container_traits< Container >
    >
class stepper_euler
{
    // provide basic typedefs
 public:

    typedef unsigned short order_type;
    typedef Time time_type;
    typedef Traits traits_type;
    typedef typename traits_type::container_type container_type;
    typedef typename traits_type::value_type value_type;


    // private members
 private:

    container_type m_dxdt;


    // public interface
 public:

    // return the order of the stepper
    order_type order_step() const { return 1; }

    // standard constructor, m_dxdt is not adjusted
    stepper_euler( void )
        {
        }

    // contructor, which adjusts m_dxdt
    stepper_euler( const container_type &x )
    {
        adjust_size( x );
    }


    // adjust the size of m_dxdt
    void adjust_size( const container_type &x )
    {
        traits_type::adjust_size( x , m_dxdt );
    }

    // performs one step with the knowledge of dxdt(t)
    template< class DynamicalSystem >
    void do_step( DynamicalSystem &system ,
                  container_type &x ,
                  const container_type &dxdt ,
                  time_type t ,
                  time_type dt )
    {
        detail::it_algebra::increment( traits_type::begin(x) ,
                                       traits_type::end(x) ,
                                       traits_type::begin(dxdt) , 
                                       dt );
    }

    // performs one step
    template< class DynamicalSystem >
    void do_step( DynamicalSystem &system ,
                  container_type &x ,
                  time_type t ,
                  time_type dt )
    {
        system( x , m_dxdt , t );
        do_step( system , x , m_dxdt , t , dt );
    }
};











template< class Container >
struct container_traits
{
    typedef Container container_type;
    typedef typename container_type::value_type value_type;
    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

    static void resize( const container_type &x , container_type &dxdt ) { dxdt.resize( x.size() ); }
    static bool same_size( const container_type &x1 , const container_type &x2 ) { return (x1.size() == x2.size()); }
    static void adjust_size( const container_type &x1 , container_type &x2 ) { if( !same_size( x1 , x2 ) ) resize( x1 , x2 ); }

    static iterator begin( container_type &x ) { return x.begin(); }
    static const_iterator begin( const container_type &x ) { return x.begin(); }
    static iterator end( container_type &x ) { return x.end(); }
    static const_iterator end( const container_type &x ) { return x.end(); }
};

