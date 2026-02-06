# Kalman_Filter_Code
## Overview
This project implements a discrete-time Kalman Filter from scratch in C++.

## Goals
- Understand Kalman filtering mathematically
- Avoid black-box linear algebra libraries
- Study numerical stability and covariance behavior

## State Definition
x = [x, y, z]^T

## Algorithm
1. Simulate true state (x_true)
2. Predict estimate (x̂⁻)
3. Update with noisy measurement
4. Track estimation error and covariance trace

## Results
- Stable estimation
- Decreasing covariance trace
- Bounded stochastic error
