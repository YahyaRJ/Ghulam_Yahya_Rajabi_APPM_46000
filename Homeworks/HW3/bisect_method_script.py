################################################################################
# This python script presents examples regarding the bisection method for
# 1D nonlinear root-finding, as presented in class.
# APPM 4650 Fall 2021
################################################################################
# Import libraries
import numpy as np;
import matplotlib.pyplot as plt;

# First, we define our example function f(x)=x+cos(x)-3.
def fun(x):
    return x + np.cos(x)-3;

# We plot our function and the x-axis y=0. We notice there is a unique intersection
# between 3 and 4.
x = np.arange(3,4,0.01);
y = fun(x);
plt.plot(x,y,x,0*x);
plt.xlabel('x'); plt.ylabel('y=f(x)');
plt.show();

input(); #pauses until user input

################################################################################
# We now implement the bisection method. Note that this can also be written in a
# separate python file and imported as a module
def bisect_method(f,a,b,tol,nmax,vrb=False):
    #Bisection method applied to f between a and b

    # Initial values for interval [an,bn], midpoint xn
    an = a; bn=b; n=0;
    xn = (an+bn)/2;
    # Current guess is stored at rn[n]
    rn=np.array([xn]);
    r=xn;
    ier=0;

    if vrb:
        print("\n Bisection method with nmax=%d and tol=%1.1e\n" % (nmax, tol));

    # The code cannot work if f(a) and f(b) have the same sign.
    # In this case, the code displays an error message, outputs empty answers and exits.
    if f(a)*f(b)>=0:
        print("\n Interval is inadequate, f(a)*f(b)>=0. Try again \n")
        print("f(a)*f(b) = %1.1f \n" % f(a)*f(b));
        r = None;
        return r
    else:
        #If f(a)f(b), we proceed with the method.
        if vrb:
            print("\n|--n--|--an--|--bn--|----xn----|-|bn-an|--|---|f(xn)|---|");

            # We start two plots as subplots in the same figure.
            fig, (ax1, ax2) = plt.subplots(1, 2); #Creates figure fig and subplots
            fig.suptitle('Bisection method results'); #Sets title of the figure
            ax1.set(xlabel='x',ylabel='y=f(x)'); #Sets labels for axis for subplot 1
            # We plot y=f(x) on the left subplot.
            xl=np.linspace(a,b,100,endpoint=True); yl=f(xl);
            ax1.plot(xl,yl);

        while n<=nmax:
            if vrb:
                print("|--%d--|%1.4f|%1.4f|%1.8f|%1.8f|%1.8f|" % (n,an,bn,xn,bn-an,np.abs(f(xn))));

                ################################################################
                # Plot results of bisection on subplot 1 of 2 (horizontal).
                xint = np.array([an,bn]);
                yint=f(xint);
                ax1.plot(xint,yint,'ko',xn,f(xn),'rs');
                ################################################################

            # Bisection method step: test subintervals [an,xn] and [xn,bn]
            # If the estimate for the error (root-xn) is less than tol, exit
            if (bn-an)<2*tol: # better test than np.abs(f(xn))<tol
                ier=1;
                break;

            # If f(an)*f(xn)<0, pick left interval, update bn
            if f(an)*f(xn)<0:
                bn=xn;
            else:
                #else, pick right interval, update an
                an=xn;

            # update midpoint xn, increase n.
            n += 1;
            xn = (an+bn)/2;
            rn = np.append(rn,xn);

    # Set root estimate to xn.
    r=xn;

    if vrb:
        ########################################################################
        # subplot 2: approximate error log-log plot
        e = np.abs(r-rn[0:n]);
        #length of interval
        ln = (b-a)*np.exp2(-np.arange(0,e.size));
        #log-log plot error vs interval length
        ax2.plot(-np.log2(ln),np.log2(e),'r--');
        ax2.set(xlabel='-log2(bn-an)',ylabel='log2(error)');
        ########################################################################

    return r, rn;
################################################################################
# Now, we run this method for our example function
(r,rn)=bisect_method(fun,3,4,5e-16,100,True);
plt.show();
input(); #pause until user input.

# Finally, we can check that |r-rn| is proportional to (1/2)^n, and so,
# log10|r-rn| ~ n (log10(1/2)) + C.
en = np.abs(r-rn[0:len(rn)-1]);
c1 = np.polyfit(np.arange(0,len(en)),np.log10(en),1);
print(10**c1); # 10 to the power of this coefficient is approximately equal to 0.5
