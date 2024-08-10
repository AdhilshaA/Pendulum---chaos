$m_1$ has 2 arms each of Length - 12”, Width - 1.5”,
Thickness - 0.375”.

$m_2$ has 1 arm each of Length - 9, Width = 1.5”, Thickness - 0.375”.



1. **Mass and Length Approximations:**
   $$
   \begin{align*}
   m_1 &= 2 \times 12 \times 1.5 \times 0.375 \times 0.04440894 \quad \text{(density of aluminum} = 0.04440894 \, \text{kg/in}^3) \\
   m_2 &= 9\times 1.5\times 0.375\times 0.04440894 \\
   l_{\text{upper}} &= 12 \\
   l_{\text{lower}} &= 9 \\
   l_3 &= l_{\text{upper}} \quad \text{(same length of both arms in the upper half)} \\
   l_2 &= \frac{l_{\text{lower}}}{2} \quad \text{(length of arm till the center of mass)} \\
   l_1 &= \frac{l_{\text{upper}}}{2} \quad \text{(length of arms till the center of mass)} \\
   g &= 9.80665 \quad \text{(acceleration due to gravity)} \\
   \end{align*}
   $$

The double pendulum is much studied system that acts
as host to a nonlinear dynamical behaviour, despite its
simple construction.

The Lyapunov exponent is a quantity that can be used to study the system's sensitivity to the initial state. Mathematically, it is defined as:

$$
d(t) = e^{\lambda t} d(0)
$$

where:
- $d(t)$ represents the separation between two initially nearby trajectories at time $t$.
- $\lambda$ is the Lyapunov exponent, which characterizes the exponential rate of divergence or convergence of nearby trajectories.
- $d(0)$ is the initial separation between the trajectories at $t = 0$.


In physics, chaos theory deals with the behavior of dynamic
systems that exhibit sensitive dependence on initial
conditions, meaning that small changes in the initial
conditions can lead to vastly different outcomes.

my coupled equations
Certainly! Here are the Python functions converted into Markdown script equations with proper fraction formatting:

Certainly! Here are the rewritten equations with subscripts and brackets only for the \(\mu1\) equation:

$$ J_a = \frac{1}{3} m_1 L_1^2 + m_2 L_1^2 $$

$$ J_b = \frac{1}{3} m_2 L_2^2 $$

$$ J_c = \frac{1}{2} m_2 L_1 L_2 $$

$$ \mu_1 = \left( \frac{1}{2} m_1 + m_2 \right) g L_1 $$

$$ \mu_2 = \frac{1}{2} m_2 g L_2 $$

The coupled equations are:

   $$ \frac{d\theta_1}{dt} = \dot{\theta_1} $$


   $$ \frac{d\theta_2}{dt} = \dot{\theta_2} $$


   $$ \frac{d\dot{\theta_1} }{dt}= \frac{ -\frac{J_c^2 }{ J_b}\dot{\theta_1}^2 \cos(\theta_2 - \theta_1) \sin(\theta_2 - \theta_1) + \frac{J_c \mu_1 }{ J_b} \cos(\theta_2 - \theta_1) \sin\theta_1 - J_c \dot{\theta_2}^2\sin(\theta_2 - \theta_1) - \mu_1 \sin\theta_1}{ J_a -  \frac{J_c^2 }{ J_b} \cos^2(\theta_2 - \theta_1) } $$


   $$ \frac{d\dot{\theta_2} }{dt}= \frac{ -\frac{J_c^2 }{ J_a}\dot{\theta_2}^2 \cos(\theta_2 - \theta_1) \sin(\theta_2 - \theta_1) + \frac{J_c \mu_1 }{ J_a} \cos(\theta_2 - \theta_1) \sin\theta_1 - J_c \dot{\theta_1}^2\sin(\theta_2 - \theta_1) - \mu_2 \sin\theta_2}{ J_b -  \frac{J_c^2 }{ J_a} \cos^2(\theta_2 - \theta_1) } $$