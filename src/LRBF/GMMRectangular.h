#ifndef GMMRECTANGLE_H
#define GMMRECTANGLE_H

// include header file
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <vector>

using namespace std;

class GMMRECTANGLE {
private:
  int LengthNode, HeightNode;
  double Length, Height;
  double LengthCenter, HeightCenter;
  double Gamma;

  vector<vector<double>> BouNodes, SouNodes;

public:
  int Nbou, Nsou, Nsurf;
  GMMRECTANGLE(int mm, int nn, double LLength, double HHeight,
               double LLengthCenter = 0.0, double HHeightCenter = 0.0,
               double GGamma = 1.5);
  ~GMMRECTANGLE(){};

  void CreateBoundaryNode();
  void CreateSourceNode();

  vector<double> GetBoundaryNode(int i) const;
  vector<double> GetSourceNode(int i) const;
  vector<double> GetSurfaceNode(int i) const;
  vector<double> GetOutNormal(int i) const;
  void UpateBoundaryNode(const vector<double> EtaDiff);
  // vector<double> GetBoundaryOutNormal(int i) const;
};

GMMRECTANGLE::GMMRECTANGLE(int mm, int nn, double LLength, double HHeight,
                           double LLengthCenter, double HHeightCenter,
                           double GGamma) {
  LengthNode = mm;
  HeightNode = nn;
  Length = LLength;
  Height = HHeight;
  Gamma = GGamma;

  LengthCenter = LLengthCenter;
  HeightCenter = HHeightCenter;

  Nbou = 2 * LengthNode + 2 * HeightNode;
  Nsou = 2 * LengthNode + 2 * HeightNode;
  Nsurf = LengthNode;

  CreateBoundaryNode();
  CreateSourceNode();

  cout << "Rectangular node distribution is created" << endl;
  cout << "Nbou: " << Nbou << ", Nsou: " << Nsou << endl;
  cout << "Gamma: " << Gamma << endl;

  assert(Gamma > 1);
}

// create coordinate of boundary nodes
void GMMRECTANGLE::CreateBoundaryNode() {
  vector<double> temp(2);
  BouNodes.resize(Nbou);

  for (int i = 0; i < Nbou; i++) {
    BouNodes[i].resize(2);

    if (i < LengthNode) {
      temp[0] = (double)(i)*Length / (double)(LengthNode) +
                .5 * Length / (double)(LengthNode);
      temp[1] = Height;
    } else if (i < LengthNode + HeightNode) {
      temp[0] = Length;
      temp[1] = Height - .5 * Height / (double)(HeightNode) -
                (double)(i - LengthNode) * Height / (double)(HeightNode);
    } else if (i < LengthNode + HeightNode + LengthNode) {
      temp[0] =
          Length - .5 * Length / (double)(LengthNode) -
          (double)(i - LengthNode - HeightNode) * Length / (double)(LengthNode);
      temp[1] = 0.;
    } else {
      temp[0] = 0.;
      temp[1] = (double)(i - LengthNode - HeightNode - LengthNode) * Height /
                    (double)(HeightNode) +
                .5 * Height / (double)(HeightNode);
    }

    temp[0] = temp[0] - Length / 2. + LengthCenter;
    temp[1] = temp[1] - Height / 2. + HeightCenter;

    BouNodes[i][0] = temp[0];
    BouNodes[i][1] = temp[1];
  }
}

// create coordinate of source nodes
void GMMRECTANGLE::CreateSourceNode() {
  assert(Nsou == Nbou);

  vector<double> temp(2);
  SouNodes.resize(Nsou);

  string type = "Circular";
  if (type == "Circular") {
    double radius, theta;
    radius = Gamma / 2 * sqrt(Length * Length + Height * Height);

    for (int i = 0; i < Nsou; i++) {
      SouNodes[i].resize(2);
      theta = 2 * M_PI * ((double)i / (double)Nsou);
      SouNodes[i][0] = radius * cos(theta) + LengthCenter;
      SouNodes[i][1] = radius * sin(theta) + HeightCenter;
    }
  } else if (type == "DomainDependent") {
    for (int i = 0; i < Nsou; i++) {
      SouNodes[i].resize(2);

      SouNodes[i][0] = (BouNodes[i][0] - LengthCenter) * Gamma + LengthCenter;
      SouNodes[i][1] = (BouNodes[i][1] - HeightCenter) * Gamma + HeightCenter;
    }
  }
}

// get coordinate of boundary nodes (index start from 1)
vector<double> GMMRECTANGLE::GetBoundaryNode(int i) const {
  assert(i >= 1 && i <= Nbou);

  return BouNodes[i - 1];
}

// get coordinate of source nodes (index start from 1)
vector<double> GMMRECTANGLE::GetSourceNode(int i) const {
  assert(i >= 1 && i <= Nsou);

  return SouNodes[i - 1];
}

// get coordinate of source nodes (index start from 1)
vector<double> GMMRECTANGLE::GetSurfaceNode(int i) const {
  assert(i >= 1 && i <= LengthNode);

  return BouNodes[i - 1];
}

// get unit out normal vectors (index start from 1, only retangular shape)
vector<double> GMMRECTANGLE::GetOutNormal(int i) const {
  assert(i >= 1 && i <= Nbou);

  vector<double> temp(2);

  if (i <= LengthNode) {
    temp[0] = 0.;
    temp[1] = 1.;
  } else if (i <= LengthNode + HeightNode) {
    temp[0] = 1.;
    temp[1] = 0.;
  } else if (i <= 2 * LengthNode + HeightNode) {
    temp[0] = 0.;
    temp[1] = -1.;
  } else {
    temp[0] = -1.;
    temp[1] = 0.;
  }

  return temp;
}

void GMMRECTANGLE::UpateBoundaryNode(vector<double> EtaDiff) {
  assert(gmm::vect_size(EtaDiff) == LengthNode);

  double Height;

  for (int i = 0; i < Nbou; i++) {

    if (i < LengthNode) {
      BouNodes[i][1] = BouNodes[i][1] + EtaDiff[i];
    } else if (i < LengthNode + HeightNode) {

      Height =
          BouNodes[LengthNode - 1][1] - BouNodes[LengthNode + HeightNode][1];

      BouNodes[i][1] =
          (BouNodes[LengthNode + HeightNode][1] + Height -
           .5 * Height / (double)(HeightNode) -
           (double)(i - LengthNode) * Height / (double)(HeightNode));

    } else if (i < LengthNode + HeightNode + LengthNode) {
      BouNodes[i][1] = BouNodes[i][1];
    } else {

      Height = BouNodes[0][1] - BouNodes[2 * LengthNode + HeightNode - 1][1];

      BouNodes[i][1] = (BouNodes[2 * LengthNode + HeightNode - 1][1] +
                        (double)(i - 2 * LengthNode - HeightNode) * Height /
                            (double)(HeightNode) +
                        .5 * Height / (double)(HeightNode));
    }
  }
}

#endif
