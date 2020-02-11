/* Boost definition.cpp

 Copyright 2009 Karsten Ahnert
 Copyright 2009 Mario Mulansky
 Copyright 2009 Andre Bergner

 The code show the usage of odeint with the Lorenz system.

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
*/

#include <iostream>
#include <iterator>

#include <boost/numeric/odeint.hpp>
#include <boost/lambda/lambda.hpp>

#define tab "\t"

using namespace std;
using namespace boost::numeric::odeint;
using boost::lambda::_1;
using boost::lambda::_2;
using boost::lambda::var;

typedef vector< double > state_type;

const double sigma = 10.0;
const double R = 28.0;
const double b = 8.0 / 3.0;

void lorenz( state_type &x , state_type &dxdt , double t )
{                                                         
    dxdt[0] = sigma * ( x[1] - x[0] );
    dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
    dxdt[2] = x[0]*x[1] - b * x[2];
}

class lorenz_class
{
 public:
    void operator()( state_type &x , state_type &dxdt , double t )
    {
        dxdt[0] = sigma * ( x[1] - x[0] );
        dxdt[1] = R * x[0] - x[1] - x[0] * x[2];
        dxdt[2] = x[0]*x[1] - b * x[2];
    }
};

int main( int argc , char **argv )
{
    const double dt = 0.01;

    // basic usage
    state_type x(3);
    x[0] = 1.0 ;
    x[1] = 0.0 ;
    x[2] = 0.0;
    stepper_euler< state_type > stepper;
    stepper.adjust_size( x );

    double t = 0.0;
    for( size_t oi=0 ; oi<10000 ; ++oi,t+=dt )
    {
        stepper.do_step( lorenz , x , t , dt );
        cout << x[0] << tab << x[1] << tab << x[2] << endl;
    }

    // use lorenz_class
    lorenz_class lorenzo;
    stepper.do_step( lorenzo , x , t , dt );


    // how to use the integrate function
    integrate_const( stepper , lorenz , x , 0.0 , 10.0 , dt , cout << _1 << tab << _2[0] << "\n" );

    vector< double > lx;
    back_insert_iterator< vector<double> > iter( lx );
    integrate_const( stepper , lorenz , x , 0.0 , 10.0 , dt , var(*iter++) = _2[0] );

    clog << endl << endl << lx.size() << endl;

    return 0;
}
