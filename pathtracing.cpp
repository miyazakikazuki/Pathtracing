#define EIGEN_NO_DEBUG // コード内のassertを無効化．
#define EIGEN_DONT_PARALLELIZE // 並列を無効化．
#define EIGEN_MPL2_ONLY // LGPLライセンスのコードを使わない．
#define PI 3.141592653589793;

#include <iostream>
#include <random>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
using namespace Eigen;
const int N = 10;

struct line{
  Vector3d org;
  Vector3d dir;
};

struct plane {
  Vector3d point;
  Vector3d normal;
  Vector2d area;
};

struct sphere {
  Vector3d center;
  double radius;
};

struct forCollision {
	double t;
	int flag;
};

class Material{
private:
	int num;
public:
	void addPlane(plane p);
	void addSphere(sphere s);
	forCollision calcCollison();
};

class Plane {

};



MatrixXd Pathtracing(
  MatrixXd input,
  struct plane eye,
  struct plane lens,
  struct plane light
);

double Le(Vector3d x, Vector3d wo);
double BSDF(Vector3d x, Vector3d wi, Vector3d wo);






int main(int argc, char const *argv[]) {
  /* code */
  MatrixXd input = MatrixXd::Zero(100,100);
  struct plane eye, lens, light;
  eye.point << 5.0, 0.0, 0.0;
  eye.area << 0.1, 0.1;
  eye.normal = -eye.point.normalized();
  lens.point << 3.2, 0.0, 0.0;
  lens.area << 1.0, 1.0;
  lens.normal = -lens.point.normalized();
  light.point << 0.0, 4.0, 3.0;
  light.area << 1.0, 1.0;
  light.normal = -light.point.normalized();
  input = Pathtracing(input, eye, lens, light);
  //std::cout << "input\n" << input << std::endl;

  double* data = new double[input.cols()*input.rows()];
  cv::Mat img(input.cols(), input.rows(), CV_64FC1, data);
  for (int i = 0; i < input.rows(); i++) {
    for (int j = 0; j < input.cols(); j++) {
      data[i * input.cols() + j] = input(input.cols() - i - 1,input.rows() - j);
    }
  }
  //std::cout << "img\n" << img << std::endl;

  cv::namedWindow("image", CV_WINDOW_AUTOSIZE|CV_WINDOW_FREERATIO);
  cv::imshow("image", img);
  cv::waitKey(0);

}

MatrixXd Pathtracing(
  MatrixXd input,
  struct plane eye,
  struct plane lens,
  struct plane light)
{
  std::random_device rnd;
  std::mt19937 mt(rnd());
  std::uniform_real_distribution<double> point(0.0, 1.0);

  Vector3d xp, x0;
  Vector3d x, tangent;
  tangent.x() = eye.normal.z() / sqrt(eye.normal.x() * eye.normal.x() + eye.normal.z() * eye.normal.z());
  tangent.y() = 0;
  tangent.z() = - eye.normal.x() / sqrt(eye.normal.x() * eye.normal.x() + eye.normal.z() * eye.normal.z());
  Vector3d binormal = eye.normal.cross(tangent);
  struct line ray = { Vector3d::Zero(), Vector3d::Zero() };
  struct line shadowray = { Vector3d::Zero(), Vector3d::Zero() };
  //オブジェクトの設定
  struct plane horizontal = { Vector3d::Zero(), Vector3d::Zero(),  Vector2d::Zero()};
  horizontal.point << 0.0, -1.0, 0.0;
  horizontal.normal << 0.0, 1.0, 0.0;
  horizontal.area << 1.0, 1.0;
  struct sphere unitsphere = { Vector3d::Zero(), 1.0 };
  double A, B, C, D;

  Vector3d raynormal, wo, wi, xl;

  double t = 0.0, tplane, tsphere, tlight,paxp,pax0,paxl;
  int paramflag;

  for (int  i = 0; i < N;  i++) {
    for (int j = 0; j < input.rows(); j++){
      for(int k = 0; k < input.cols(); k++){
        xp = eye.point + ((k + point(mt)) / input.cols() - 0.5) * eye.area.x() * tangent
              + ((j + point(mt)) / input.rows() - 0.5) * eye.area.y() * binormal;

         paxp = input.cols() * input.rows() / (eye.area.x() * eye.area.y());

        /*レンズの一点をサンプリング*/
        x0 = lens.point + ((k + point(mt)) / input.cols() - 0.5) * lens.area.x() * tangent
              + ((j + point(mt)) / input.rows() - 0.5) * lens.area.y() * binormal;
        pax0 = input.cols() * input.rows() / (lens.area.x() * lens.area.y());

        ray.org = x0;
        ray.dir = (x0 - xp).normalized();
        double psigma = 1.0 / 2.0 * PI;
        double alpha = 1.0;//ray.dir.dot(eye.normal) / pax0 * psigma;/*重みづけ未考慮*/
        //std::cout << "start" << std::endl;


        do{//ray をトレースして衝突がある
		  t = 0;
          tplane = (horizontal.point - ray.org).dot(horizontal.normal) / ray.dir.dot(horizontal.normal);
          //if ((ray.org + tplane * ray.dir - horizontal.point).norm() > 5.0) {
          //  tplane = -1.0;
          //}
          tlight = (light.point - ray.org).dot(light.normal) / ray.dir.dot(light.normal);
          if ((ray.org + tlight * ray.dir - light.point).norm() > 1.0) {
            tlight = -1.0;
          }

          A = pow(ray.dir.norm(), 2);
          B = ray.dir.dot(ray.org - unitsphere.center);
          C = pow((ray.org - unitsphere.center).norm(), 2) - pow(unitsphere.radius, 2);
          D = B * B - A * C;

          if (D < 0) {
            tsphere = -1.0;
          }
          else {
            tsphere = (-B - sqrt(D)) / A;
          }
          //std::cout << tsphere << std::endl;

          tplane = 1.0 / tplane;
          tlight = 1.0 / tlight;
          tsphere = 1.0 / tsphere;

          if (tlight < tplane) {
            if (tplane < tsphere) {
              t = tsphere;
              paramflag = 2;
            }
            else {
              t = tplane;
              paramflag = 0;
            }
          }
          else if (tlight < tsphere) {
            t = tsphere;
            paramflag = 2;
          }
          else {
            t = tlight;
            paramflag = 3;
          }
          t = fmax(t, 0);

          //std::cout <<  t << std::endl;
          if (t <= 0) break;

          x = ray.org + ray.dir / t ;
          wo = -ray.dir;
          switch (paramflag) {
            case 0:
              raynormal = horizontal.normal;
              break;
            case 1:
            case 2:
              raynormal = x - unitsphere.center;
              break;
            case 3:
              input(j, k) = input(j, k) + alpha * Le(x, wo);
              //std::cout << input(j, k) << std::endl;
              break;
          }
          if (paramflag == 3) {
            //std::cout << "end1" << std::endl;
            break;
          }
          wi << point(mt), point(mt), point(mt);
          if (wi.dot(raynormal) < 0) {
            wi = -wi;
          }
          wi = wi / wi.norm();

          double psigmawi = 1.0;

		  xl = light.point + ((k + point(mt)) / input.cols() - 0.5) * light.area.x() * tangent
			  + ((j + point(mt)) / input.rows() - 0.5) * light.area.y() * binormal;

		  paxl = 1.0;

		  shadowray.org = x;

		  shadowray.dir = (xl - x).normalized();

		  if ((xl - x).dot(raynormal) > 0) {

			  input(j, k) = input(j, k) + alpha * Le(xl, x) * BSDF(x, wi, wo) / paxl;

		  }

          alpha = alpha * BSDF(x, wi, wo) * wi.dot(raynormal) / psigmawi;

          double prr = 0.5;

          if(point(mt) > prr){
            //std::cout << "end2" << std::endl;
            break;
          }
          ray.org = x;
          ray.dir = wi;

          alpha = alpha / prr;
       } while (t > 0.0);
      }
    }
    std::cout << "loading..." << (double)(i + 1) / (double)N * 100.0 << "%" << std::endl;
  }

  return input / (double)N;
}

double Le(Vector3d x, Vector3d wo){
  return 1.0;
}

double BSDF(Vector3d x, Vector3d wi, Vector3d wo){
  return 1.0;
}
