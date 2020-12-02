template <int dim>
  class InitialValuesU : public Function<dim>
  {
  public:
    InitialValuesU () : Function<dim>() {}

    virtual double value (const Point<dim>   &/*p*/,
                          const unsigned int  component = 0) const;
  };


  template <int dim>
  class InitialValuesV : public Function<dim>
  {
  public:
    InitialValuesV () : Function<dim>() {}

    virtual double value (const Point<dim>   &/*p*/,
                          const unsigned int  component = 0) const;
  };



	template <int dim>
  class InitialValueswave_speed : public Function<dim>
  {
  public:
    InitialValueswave_speed () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };

	template <int dim>
  double InitialValueswave_speed<dim>::value (const Point<dim>  &p,
                                     const unsigned int component) const
  {
    (void) component;
    Assert(component == 0, ExcIndexRange(component, 0, 1));
    double   radius=0.25;
		  if (p.square() < radius*radius)
			 return 0;
		  else
			 return 1;
  }








  template <int dim>
  double InitialValuesU<dim>::value (const Point<dim>  &/*p*/,
                                     const unsigned int component) const
  {
    (void) component;
    Assert(component == 0, ExcIndexRange(component, 0, 1));
    return 0;
  }



  template <int dim>
  double InitialValuesV<dim>::value (const Point<dim>  &/*p*/,
                                     const unsigned int component) const
  {
    (void) component;
    Assert(component == 0, ExcIndexRange(component, 0, 1));
    return 0;
  }

	









  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };






  template <int dim>
  double RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                                    const unsigned int component) const
  {
    (void) component;
    Assert(component == 0, ExcIndexRange(component, 0, 1));
    return 0;
  }



  template <int dim>
  class BoundaryValuesU : public Function<dim>
  {
 
  public:
    BoundaryValuesU () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };




  template <int dim>
  class BoundaryValuesV : public Function<dim>
  {
  public:
    BoundaryValuesV () : Function<dim>() {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
  };




  template <int dim>
  double BoundaryValuesU<dim>::value (const Point<dim> & p,
                                      const unsigned int component) const
  {
    (void) component;
    Assert(component == 0, ExcIndexRange(component, 0, 1));
	 double N=10.0;
	 
     if ((this->get_time() <= 1.0/N) && (p[0] ==-1.0) ) 
        return std::sin(this->get_time() * N * numbers::PI);
     else
      return 0;
  }


 

  template <int dim>
  double BoundaryValuesV<dim>::value (const Point<dim> & p,
                                      const unsigned int component) const
  {
    (void) component;
    Assert(component == 0, ExcIndexRange(component, 0, 1));
		double N=10.0;			
		if ((this->get_time() <= 1.0/N) && (p[0] ==-1.0) ) 
         return (std::cos(this->get_time() * N * numbers::PI) * N * numbers::PI);
      
      else
      return 0;
  }

template <int dim>
double rho(const Point<dim> &p)
{
 
  double   radius=0.25;
  if (p.square() < radius*radius)
    return 1;
  else
    return 1;
}

template <int dim>
double a(const Point<dim> &p)
{

   
  double  radius=std::sqrt(p.square()),speed,c1=0.,c2=1.,epsilon=10.0*.01,a=0.25;
  
  speed=c1+(c2-c1)/(1+std::exp(-(radius-a)/epsilon));
  
    return speed;
}
