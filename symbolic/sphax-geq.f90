geqrhs(0) = (-(2*vel(2)*vel(3)*(pot(2)*dpot(5,2)-dpot(2,2)*pot(5)&
&))/(2*pot(1)*pot(5)+2*pot(2)**2))-(2*vel(1)*vel(3)*(pot(2)*dpot(5&
&,1)-dpot(2,1)*pot(5)))/(2*pot(1)*pot(5)+2*pot(2)**2)-(2*vel(0)*ve&
&l(2)*(dpot(1,2)*pot(5)+pot(2)*dpot(2,2)))/(2*pot(1)*pot(5)+2*pot(&
&2)**2)-(2*vel(0)*vel(1)*(dpot(1,1)*pot(5)+pot(2)*dpot(2,1)))/(2*p&
&ot(1)*pot(5)+2*pot(2)**2)
 
geqrhs(1) = ((vel(3)**2*dpot(5,1))/pot(3))/2.d+0+((vel(2)**2*dpot&
&(4,1))/pot(3))/2.d+0-(vel(1)*vel(2)*dpot(3,2))/pot(3)-((vel(1)**2&
&*dpot(3,1))/pot(3))/2.d+0+(vel(0)*dpot(2,1)*vel(3))/pot(3)-((vel(&
&0)**2*dpot(1,1))/pot(3))/2.d+0
 
geqrhs(2) = ((vel(3)**2*dpot(5,2))/pot(4))/2.d+0-((vel(2)**2*dpot&
&(4,2))/pot(4))/2.d+0-(vel(1)*vel(2)*dpot(4,1))/pot(4)+((vel(1)**2&
&*dpot(3,2))/pot(4))/2.d+0+(vel(0)*dpot(2,2)*vel(3))/pot(4)-((vel(&
&0)**2*dpot(1,2))/pot(4))/2.d+0
 
geqrhs(3) = (-(2*vel(2)*vel(3)*(pot(1)*dpot(5,2)+pot(2)*dpot(2,2)&
&))/(2*pot(1)*pot(5)+2*pot(2)**2))-(2*vel(1)*vel(3)*(pot(1)*dpot(5&
&,1)+pot(2)*dpot(2,1)))/(2*pot(1)*pot(5)+2*pot(2)**2)-(2*vel(0)*ve&
&l(2)*(pot(1)*dpot(2,2)-dpot(1,2)*pot(2)))/(2*pot(1)*pot(5)+2*pot(&
&2)**2)-(2*vel(0)*vel(1)*(pot(1)*dpot(2,1)-dpot(1,1)*pot(2)))/(2*p&
&ot(1)*pot(5)+2*pot(2)**2)
 
