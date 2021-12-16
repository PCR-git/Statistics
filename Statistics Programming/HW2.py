# -*- coding: utf-8 -*-
"""

 Stat 202A 2019 Fall - Homework 02
 Author: Peter Racioppo
 Date : 10/13/2019

 INSTRUCTIONS: Please fill in the corresponding function. Do not change function names, 
 function inputs or outputs. Do not write anything outside the function.
 
"""

import numpy as np

### Part 1 : Metropolis Algorithm

def sample_uniform(size=10000, low=0, high=1):
    """
    Function 1A: Uniform (Same as HW1)
    Detail : This function return a np.ndarray with given size.
    """
    # Constants
    a = np.power(7,5)
    b = 0
    M = np.power(2,31)-1
    Y = np.zeros(size+5) # Initialize array, accounting for burn in of 5 points
    # Set a seed as the last 4 digits of the current time in microseconds
    import time
    # Y[0] = int(time.time()*10e6)%1000 # Set first value of array to seed
    # Trying to be very sure that the seed is robust:
    Y[0] = int((time.time()*10e6)**2)%1000 + int((time.time()*10e6))%1000
    # Generate uniform random samples:
    for i in (np.arange(size+4)+1): # (account for burn in)
        Y[i] = (a*Y[i-1] + b)%M
    U = Y[5:]/M  # Uniform over [0,1] (include burn in)
    X = U*(high-low) + low  # Uniform over [low,high]
    return X

def sample_normal_chain(x0, c, chain_length=100, mean=0, var=1):

    """
    Function 1B: Normal by Metropolis
    Detail :    This function returns multiple chains by Metropolis sampling from N(mean, var).
                For every train, the proposal distribution at x is y ~ Uniform[x-c, x+c].
                The input x0 is a 1 dimension np.ndarray. 
                Return a np.ndarray X with size x0.shape[0] * chain_length. In X, each row X[i]
                should be a chain and t-th column X[:, t] corresponding to all different X_t. 
    """
    
    X = np.zeros([x0.shape[0],chain_length]) # Initialize X
    for i in np.arange(len(x0)):
        x_vec = [x0[i]] # Starting value of MCMC
        # Sample from a uniform distribution:
        U = sample_uniform(chain_length,0,1)
        # Sample from a different uniform distribution:
        samp = sample_uniform(chain_length,0,1)
        for j in np.arange(chain_length-1):
            x = x_vec[-1] # X_t
            y = samp[j]*2*c + (x-c) # X_t+1 (sample from proposal distrb.)
            r = np.exp(((x-mean)**2-(y-mean)**2)/(2*var)) # p(X_t+1)/p(X_t)
            if np.minimum(1, r) > U[j]:
                x = y
            x_vec.append(x)
        X[i,:] = x_vec

    return X

def metropolis_simulation(num_chain=1000, chain_length=100, mean=0, var=1):
    """
    Function 1C: Simulate metropolis with different setting.
    Detail: Try different settings and output movies of histograms.
    """
    import math
    list_a = [0, -0.1, -1, -10] # Lower limit for starting sample
    list_b = [1, 0.1, 1, 10] # Upper limit for starting sample
    list_c = [1, 2] # Sampling parameter
    # list_a = [-1]
    # list_b = [1]
    # list_c = [1]

    for a, b, c in [(a, b, c) for a in list_a for b in list_b for c in list_c]:
        # --------------------
        # Metropolis Algorithm

        # Sample the starting point from a uniform distribution:
        x0 = sample_uniform(size=num_chain,low=a,high=b)

        # Run the Metropolis Algorithm:
        X_normal = sample_normal_chain(x0,c,chain_length,mean,var)

        burn_in = math.ceil(0.5*chain_length) # Burn in
        X_f = X_normal[:,burn_in:].flatten() # Remove values for burn in

        # Plot movie and save 
        # Here plot chain_length graphs, each of them is a histogram of num_chain point. 
        # You may use matplotlib.animation and matplotlib.rc to save graphs into gif movies.

        # --------------------
        # Histogram

        # Import plot stuff:
        import matplotlib.pyplot as plt
        import matplotlib
        from pylab import rcParams
        matplotlib.rcParams.update({'font.size':16})

        # Histogram of the samples:
        sigma = np.sqrt(var)
        fig, ax = plt.subplots()
        lim = math.ceil(max(abs(np.array(X_f))))
        num_bins = math.ceil((lim/6)*100)
        n, bins, patches = plt.hist(X_f, num_bins,edgecolor='black',density='normed',alpha=0.5,linewidth=2,color='b')

        # Overlaying the target PDF:
        import scipy.stats as stats
        x = np.linspace(mean - lim*sigma, mean + lim*sigma, 500)
        y = stats.norm.pdf(x, mean, sigma)
        plt.plot(x,y,color='black',linestyle='--',linewidth=2)
        # Plot stats
        str1 = 'Histogram of Metropolis Samples: [a,b,c]='
        str2 = str([a,b,c])
        title = str1+str2
        fig.suptitle(title, fontsize=20)
        plt.xlabel('x', fontsize=20)
        plt.ylabel('Density', fontsize=20)
        # ax.set_xlim([-lim,lim])
        ax.set_xlim([-4,4])
        # rcParams['figure.figsize'] = 12, 6
        ax.grid()
        plt.show()
        # fig.savefig(title+'.png')

        # --------------------
        # Animation

        import matplotlib.patches as patches
        import matplotlib.path as path
        import matplotlib.animation as animation

        mean = 0 # Mean of the normal distribution
        var = 1 # Variance of the normal distribution
        sigma = np.sqrt(var)
        lim = 4
        num_bins = 60
        samp_length = chain_length - burn_in

        frac = 1 # We record 1/frac frames
        num_frames = int(len(X_f)/(samp_length*frac))

        fig = plt.figure(figsize=(10,6))

        def update_hist(num, X_f):
            plt.cla()
            n, bins, patches = plt.hist(X_f[0:num*samp_length*frac],num_bins,edgecolor='black',density='normed',\
            							alpha=0.5,linewidth=2,color='b')

            x = np.linspace(mean - lim*sigma, mean + lim*sigma, 500)
            y = stats.norm.pdf(x, mean, sigma)
            plt.plot(x,y,color='black',linestyle='--',linewidth=2)

            plt.xlim(-lim, lim)
            plt.ylim(0, 0.5)
            plt.xlabel('x',fontsize=20)
            plt.ylabel('Density',fontsize=20)
            str1 = 'Histogram of Metropolis Samples: [a,b,c]='
            str2 = str([a,b,c])
            title = str1+str2
            fig.suptitle(title, fontsize=20)
            plt.grid()

        ani = animation.FuncAnimation(fig,update_hist,num_frames,fargs=(X_f, ),repeat='False', interval=1)
        plt.show()
        # ani.save('Metropolis_Movie.mp4', fps=80)


### Part 2 : Gibbs Sampling

def gibbs_sample(x0, y0, rho, num_chain=1000, chain_length=100, mean=0, var=1):
    """
    Function 2A: Bivariate normal with correlation rho
    Detail :    This function return multiple chains by Gibbs sampling
                The input x0, y0, rho, num_chain is a number. This time, we use same starting point. 
                Return a np.ndarray X with size num_chain * chain_length * 2. In X, each row X[i]
                should be a chain and t-th column X[:, t] corresponding to all different pair (X_t, Y_t). 
    """

    def sample_normal(mean=0, var=1):
        """
        This function samples from a normal distribution.
        """
        # Draw from 2 uniform distributions
        U1 = sample_uniform(2,0,1)[0]
        U2 = sample_uniform(2,0,1)[1]
        #     import random
        #     U1 = random.uniform(0,1)
        #     U2 = random.uniform(0,1)

        # Polar coordinates
        pi = 3.14159265
        theta = 2*pi*U1
        R = np.sqrt(-2*np.log(1-U2))

        # Cartesian coordinates
        z = R*np.cos(theta)
        x = z*np.sqrt(var) + mean

        return x

    # Gibbs Sampling
    X = np.zeros([num_chain,chain_length,2]) # Initialize X
    for i in np.arange(num_chain):
        # Starting values of MCMCs:
#         x_vec = [x0[i]]
#         y_vec = [y0[i]]
        x_vec = [x0]
        y_vec = [y0]
        for j in np.arange(chain_length-1):
            x = x_vec[-1]
            y = y_vec[-1]
    #         x = np.random.normal(rho*y,1-rho**2)
    #         y = np.random.normal(rho*x,1-rho**2)
            x = sample_normal(rho*y,1-rho**2)
            y = sample_normal(rho*x,1-rho**2)
            x_vec.append(x)
            y_vec.append(y)
        X[i,:,0] = x_vec
        X[i,:,1] = y_vec
    
    return X

def gibbs_simulation():
    """
    Function 2B: Simulate Gibbs with different rho and plot
    Detail : Try different setting and output movie of histgrams. 
             Discard first 50 steps and output 50~100 steps only.
    """
    
    # Imports:
    import math
    import matplotlib.pyplot as plt
    import matplotlib
    from pylab import rcParams
    matplotlib.rcParams.update({'font.size':16})
    
    # Parameters
    mean = 0
    var = 1
    num_chain=1000
    chain_length=100
    x0 = 0
    y0 = 0
    
#     # Sample the starting point from a uniform distribution:
#     a = -1 # Lower limit for starting sample
#     b = 1 # Upper limit for starting sample
#     x0 = sample_uniform(size=1,low=a,high=b)[0]
#     y0 = sample_uniform(size=1,low=a,high=b)[0]

    list_rho = [0, 0.2, -1, 1]
    # list_rho = [0]
    for rho in list_rho:
        # Run Gibbs Sampling
        X = gibbs_sample(x0,y0,rho,num_chain,chain_length,mean,var)

        # Discard first 50 steps
        burn_in = 50 # Burn in
        X_f = X[:,burn_in:,:] # Remove values
        X_data = X_f[:,:,0].flatten()
        Y_data = X_f[:,:,1].flatten()
        
        # Plot info
        sigma = np.sqrt(var)
        fig, ax = plt.subplots()
        # lim = math.ceil(max(abs(np.array(X_data))))
        lim = 4
#         num_bins = math.ceil((lim/6)*100)
        num_bins = 80
        
        # Histogram of the X samples:
        n1, bins1, patches2 = plt.hist(X_data,num_bins,density='normed',linewidth=2,\
                                    edgecolor='black',alpha = 0.5,color='b')
        # Histogram of the Y samples:
        n2, bins2, patches2 = plt.hist(Y_data,num_bins,density='normed',linewidth=2,\
                                    edgecolor='black',alpha = 0.5,color='r')
        
        # Overlaying the target PDF:
        import scipy.stats as stats
        x = np.linspace(mean - 4*sigma, mean + 4*sigma, 500)
        y = stats.norm.pdf(x, mean, sigma)
        plt.plot(x,y,color='black',linestyle='--',linewidth=2)
        
        # Plot info
        str1 = 'Histogram of Gibbs Samples: rho='
        str2 = str(rho)
        title = str1+str2
        fig.suptitle(title, fontsize=20)
        plt.xlabel('X', fontsize=20)
        plt.ylabel('Density', fontsize=20)
        # ax.set_xlim([-lim,lim])
        ax.set_xlim([-4,4])
        rcParams['figure.figsize'] = 12, 6
        ax.legend(['Target','X','Y'])
        ax.grid()
        plt.show()
        fig.savefig(title+'.png')

        # --------------------
        if np.abs(rho)!=1:
            # Scatter plot of a single chain
            X1 = X_data[0:chain_length-burn_in]
            Y1 = Y_data[0:chain_length-burn_in]

            fig, ax = plt.subplots()
            for i in np.arange(chain_length-burn_in-1):
                plt.arrow(X1[i],Y1[i],X1[i+1]-X1[i],Y1[i+1]-Y1[i],color='b',shape='full',\
                          length_includes_head=True,head_length=0.1,head_width=0.1,linestyle='--')
                ax.scatter(X1,Y1,c='k')
            
            # Plot info
            str1 = 'Scatterplot of Single Chain, Gibbs Normal Samples: rho='
            str2 = str(rho)
            title = str1+str2
            fig.suptitle(title, fontsize=20)
            plt.xlabel('X', fontsize=20)
            plt.ylabel('Y', fontsize=20)
            # ax.set_xlim([-lim,lim])
            ax.set_xlim([-4,4])
            rcParams['figure.figsize'] = 12, 6
            ax.grid()
            plt.show()
            fig.savefig(title+'.png')

        # --------------------
        # Plot movie and save

        import matplotlib.patches as patches
        import matplotlib.path as path
        import matplotlib.animation as animation

        samp_length = chain_length - burn_in
        frac = 1 # We record 1/frac frames
        num_frames = int(len(X_data)/(samp_length*frac))

        fig = plt.figure(figsize=(12,6))

        def update_hist(num, X_data, Y_data):
            plt.cla()
            n1, bins1, patches1 = plt.hist(X_data[0:num*samp_length*frac],num_bins,density='normed',linewidth=2,\
                                    edgecolor='black',alpha = 0.5,color= 'b')
            n2, bins2, patches2 = plt.hist(Y_data[0:num*samp_length*frac],num_bins,density='normed',linewidth=2,\
                                    edgecolor='black',alpha = 0.5,color= 'r')
            
            x = np.linspace(mean - 4*sigma, mean + 4*sigma, 500)
            y = stats.norm.pdf(x, mean, sigma)
            plt.plot(x,y,color='black',linestyle='--',linewidth=2)

            plt.xlim(-lim, lim)
            plt.ylim(0, 0.5)
            plt.xlabel('X, Y',fontsize=20)
            plt.ylabel('Density',fontsize=20)
            str1 = 'Gibbs Sampling of Normal Distribution: rho='
            str2 = str(rho)
            title = str1+str2
            fig.suptitle(title, fontsize=20)
            ax.legend(['Target','X','Y'])
            plt.grid()

        ani = animation.FuncAnimation(fig,update_hist,num_frames,fargs=(X_data,Y_data),repeat='False',interval=1)
        plt.show()
        # ani.save('Gibbs_Movie.mp4', fps=80)

# TESTS:
# metropolis_simulation(num_chain=1000, chain_length=100, mean=0, var=1)
# gibbs_simulation()


