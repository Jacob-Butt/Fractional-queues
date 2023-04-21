# G/G/1 queue using Mittag-Leffler arrivals

# See https://www.r-bloggers.com/2010/05/simulating-a-queue-in-r/

# Parameters

t.end   <- 10^6 # duration of sim
t.clock <- 0    # sim time
Ta <- 1    # interarrival period
beta1 <- 0.7000	# interarrival ML parameter (1 is exponential)
Ts <- 1    # service period
beta2 <- 0.7000 # interservice ML parameter (1 is exponential)
t1 <- 0         # time for next arrival (the first event is an arrival)
t2 <- t.end     # time for next departure (the first event is an arrival)
tn <- t.clock   # tmp var for last event time
tb <- 0         # tmp var for last busy-time start
n <- 0          # number in system
s <- 0          # cumulative number-time product
b <- 0          # total busy time
c <- 0          # total completions
qc <- 0         # plot instantaneous q size
tc <- 0         # plot time delta
plotSamples <- 10000
set.seed(1) #Use this for debug purposes

# Loop

while (t.clock < t.end) {
    if (t1 < t2) {      # arrival event
        t.clock <- t1
        s <- s + n * (t.clock - tn)  # delta time-weighted number in queue
        n <- n + 1
        if (t.clock < plotSamples) { 
            qc <- append(qc,n)
            tc <- append(tc,t.clock) 
        }
        tn <- t.clock
        t1 <- t.clock -Ta*log(runif(1))*(sin(beta1*pi)/tan(beta1*pi*runif(1))-cos(beta1*pi))^(1/beta1)
        if(n == 1) { 
            tb <- t.clock
            t2 <- t.clock -Ts*log(runif(1))*(sin(beta2*pi)/tan(beta2*pi*runif(1))-cos(beta2*pi))^(1/beta2)  # interarrival period
        }
    } else {            # departure event
        t.clock <- t2
        s <- s + n * (t.clock - tn)  # delta time-weighted number in queue
        n <- n - 1
        if (t.clock < plotSamples) { 
            qc <- append(qc,n)
            tc <- append(tc,t.clock)
        }
        tn <- t.clock
        c <- c + 1
        if (n > 0) { 
            t2 <- t.clock -Ts*log(runif(1))*(sin(beta2*pi)/tan(beta2*pi*runif(1))-cos(beta2*pi))^(1/beta2)  # service period
        }
        else { 
            t2 <- t.end
            b <- b + t.clock - tb
        }
    }   
}

# Some of these variables may not make sense for G/G/1 queues. You can get M/M/1 by setting beta1=1 and beta2=1

u <- b/t.clock       # utilization B/T
N <- s/t.clock       # mean queue length (see the Load Average notes)
x <- c/t.clock       # mean throughput C/T
r <- N/x             # mean residence time (from Little's law: Q = XR)
q <- sum(qc)/max(tc) # estimated queue length for plot

#par(mar=c(5,6,4,1)+.1) # used to fit plot into the viewing window if necessary

plot(tc,qc,type="s",xlab="Time",ylab="Instantaneous queue length")


