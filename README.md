# Computer-Vision
Code and project work for Computer Vision

## Project 1 - Motion Detection via Spatial and Temporal Filtering
This project uses various algorithms which utilize
temporal derivatives in order to detect moving people using a video recorded on a stationary
camera. The zip file is provided which contains all the images that combined make a video. 
Using a one dimensional, temporal derivative filtering mask, we can
detect drastic shifts in the grayscale intensity of pixels on a stationary background
when a moving object passes them along the time spectrum. By looking at the
large temporal gradients in the pixels where a moving object passed, we can detect
and track the motion of the moving objects. To accomplish this, we created a
binary mask by thresholding the absolute value of the temporal derivatives and
combined the result with the mask to highlight the pixels around the moving
objects. Several 2D spatial smoothing algorithms were also used including two
box filters and a gaussian filter with varying the sigma value.


## Project 2 - Image Mosaicing
This projects utilizes various scripts and function to create a projective transformation
between two image spaces, otherwise known as image mosaicing. Two images are warped into the same image space 
to create an image that is the union of both images. This is accomplished by detecting corner points in a image,
correlating those points with corresponding points in a similar image, using those
points to estimate a homography between the two image spaces, and warping the
images into the same image spaces. A Harris corner detection function in MATLAB 
was implemented and used to apply normalized cross correlation to
find corresponding corner features between the given images. RANSAC was also used to estimate the homography
in a fashion robust against incorrect correspondences. Finally, images are warped together using reverse image warping.

## Project 3 - Optical Flow
In computer vision, there is sometimes a need for understanding the apparent
displacement of an object of interest between several images of a scene, whether
it is due to motion of the objects in the scene or the motion of the camera. This
information is known as the optical flow, and it is useful in many video based
tracking applications. In this project, an algorithm for estimating
the optical flow, the Lucas-Kanade algorithm, is implemented. This algorithm 
can be used to estimate the flow vectors in a benchmark set of images.

## Project 4 - Target Tracking
Target tracking is an important topic in computer vision, having application in
automatic surveillance applications and, in general, in video processing applications.
One challenge in target tracking is detecting and recovery from occlusion.
In this project, various algorithms for occlusion detection and recovery are implemented. Peak-to-sidelobe ratio is used 
as a detection statistic for occlusion detection and two methods for occlusion recovery are experimented, which both predict the
location of the occluded object based on amodel for the dynamics of the object.
The first method is a simple constant velocity model, and the other is an adaptive
model that automatically estimates the complexity of the motion with the Hankel
matrix and predicts the position of the occluded object via a linear predictor.
