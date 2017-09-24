
// 1. マクロを定義．
#define EIGEN_NO_DEBUG // コード内のassertを無効化．
#define EIGEN_DONT_PARALLELIZE // 並列を無効化．
#define EIGEN_MPL2_ONLY // LGPLライセンスのコードを使わない．
#define PI 3.141592653589793;

// 2. Eigenをインクルード．
#include <iostream>
#include <random>
#include <opencv2/core/core.hpp>
#include <opencv2/opencv.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
using namespace Eigen;
const int N = 100;

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
  eye.point << 0.0, 0.0, 5.0;
  eye.area << 0.1, 0.1;
  eye.normal << 0.0, 0.0, -1.0;
  lens.point << 0.0, 0.0, 3.0;
  lens.area << 1.0, 1.0;
  lens.normal << 0.0, 0.0, -1.0;
  light.point << 3.0, 4.0, 0.0;
  light.area << 1.0, 1.0;
  light.normal << -0.6, -0.8, 0.0;
  input = Pathtracing(input, eye, lens, light);
  //std::cout << "input\n" << input << std::endl;



  double* data = new double[input.cols()*input.rows()];
  cv::Mat img(input.cols(), input.cols(), CV_64FC1, data);
  for (int i = 0; i < input.rows(); i++) {
    for (int j = 0; j < input.cols(); j++) {
      data[i * input.cols() + j] =input(i,j);
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

  Vector3d xp = eye.point;
  Vector3d x0 = lens.point;
  Vector3d x, tangent;
  tangent.x() = - eye.normal.z() / sqrt(eye.normal.x() * eye.normal.x() + eye.normal.z() * eye.normal.z());
  tangent.y() = 0;
  tangent.z() = eye.normal.x() / sqrt(eye.normal.x() * eye.normal.x() + eye.normal.z() * eye.normal.z());
  Vector3d binormal = tangent.cross(eye.normal);
  struct line ray = { Vector3d::Zero(), Vector3d::Zero() };
  //オブジェクトの設定
  struct plane horizontal = { Vector3d::Zero(), Vector3d::Zero(),  Vector2d::Zero()};
  horizontal.point << 0.0, -1.0, 0.0;
  horizontal.normal << 0.0, 1.0, 0.0;
  horizontal.area << 1.0, 1.0;
  struct sphere unitsphere = { Vector3d::Zero(), 1.0 };

  Vector3d raynormal, wo, wi, normal;

  double t = 0, tplane, tsphere1, tsphere2, tlight;
  int paramflag;

  for (int  i = 0; i < N;  i++) {
    for (int j = 0; j < input.rows(); j++){
      for(int k = 0; k < input.cols(); k++){
        xp = eye.point + ((j + point(mt)) / input.cols() - 0.5) * eye.area.x() * tangent
              + ((k + point(mt)) / input.rows() - 0.5) * eye.area.y() * binormal;

        double paxp = 1 / (eye.area.x() * eye.area.y());

        /*レンズの一点をサンプリング*/
        x0 = lens.point + ((j + point(mt)) / input.cols() - 0.5) * lens.area.x() * tangent
              + ((k + point(mt)) / input.rows() - 0.5) * lens.area.y() * binormal;
        double pax0 = 1.0;

        ray.org = x0;
        ray.dir = (x0 - xp).normalized();
        double pasigma = 1.0 / 2.0 * PI;
        double alpha = ray.dir.dot(eye.normal) / pax0 * pasigma;/*重みづけ未考慮*/

        tplane = (horizontal.point- ray.org).dot(horizontal.normal) / ray.dir.dot(horizontal.normal);
        tlight = (light.point - ray.org).dot(light.normal) / ray.dir.dot(light.normal);
        tsphere1 = (-(ray.org - unitsphere.center).dot(ray.dir)
                    + ray.dir.norm() * sqrt(pow((ray.org - unitsphere.center).norm(), 2)  - unitsphere.radius * unitsphere.radius))
                  / ray.dir.dot(ray.dir);
        tsphere2 = (-(ray.org - unitsphere.center).dot(ray.dir)
          - ray.dir.norm() * sqrt(pow((ray.org - unitsphere.center).norm(), 2) - unitsphere.radius * unitsphere.radius))
          / ray.dir.dot(ray.dir);

        tplane = 1.0 / tplane;
        tlight = 1.0 / tlight;
        tsphere1 = 1.0 / tsphere1;
        tsphere2 = 1.0 / tsphere2;

        t = 0;

        if (tlight < tplane){
          if (tplane < tsphere1){
            if(tsphere1 < tsphere2) {
              t = tsphere2;
              normal = x;
              paramflag = 2;
            }else {
              t = tsphere1;
              normal = x;
              paramflag = 1;
            }
          }else {
            t = tplane;
            normal = horizontal.normal;
            paramflag = 0;
          }
        } else {
          paramflag = 3;
        }
        t = fmax(t, 0);


        while(t > 0){//ray をトレースして衝突がある
          x = ray.org + ray.dir / t ;
          wo = -ray.dir;
          switch (paramflag) {
            case 0:
              raynormal = horizontal.normal;
              break;
            case 1:
            case 2:
              raynormal = x;
              break;
            case 3:
              input(j, k) = input(j, k) + alpha * Le(x, wo);
              break;
          }
          if(paramflag == 3) break;
          wi << point(mt), point(mt), point(mt);
          wi = wi / wi.norm();

          double psigmawi = 1.0;

          alpha = alpha * BSDF(x, wi, wo) * wi.dot(normal) / psigmawi;

          double prr = 0.5;

          if(point(mt) > prr){
            break;
          }
          ray.org = x;
          ray.dir = wi;

          alpha = alpha / prr;


          tplane = -ray.org.dot(horizontal.normal) / ray.dir.dot(horizontal.normal);
          tlight = -ray.org.dot(light.normal) / ray.dir.dot(light.normal);
          tsphere1 = (-(ray.org - unitsphere.center).dot(ray.dir)
                      + ray.dir.norm() * sqrt(pow((ray.org - unitsphere.center).norm(), 2)  - unitsphere.radius * unitsphere.radius))
                    / ray.dir.dot(ray.dir);
          tsphere2 = (-(ray.org - unitsphere.center).dot(ray.dir)
            - ray.dir.norm() * sqrt(pow((ray.org - unitsphere.center).norm(), 2) - unitsphere.radius * unitsphere.radius))
            / ray.dir.dot(ray.dir);

          tplane = 1 / tplane;
          tlight = 1 / tlight;
          tsphere1 = 1 / tsphere1;
          tsphere2 = 1 / tsphere2;

          if (tlight < tplane){
            if (tplane < tsphere1){
              if(tsphere1 < tsphere2) {
                t = tsphere2;
                paramflag = 2;
              }else {
                t = tsphere1;
                paramflag = 1;
              }
            }else {
              t = tplane;
              paramflag = 0;
            }
          } else {
            paramflag = 3;
          }
          t = fmax(t, 0);

          //input(1,1) = 1;
        }
      }
    }
  }
  return input / (double)N;
}

double Le(Vector3d x, Vector3d wo){
  return 1.0;
}

double BSDF(Vector3d x, Vector3d wi, Vector3d wo){
  return 1.0;
}
