
#include "stats.h"

stats::stats(PNG & im){
  for (int x = 0; x < im.width(); x++) {
    vector<double> vHX;
    vector<double> vHY;
    vector<double> vS;
    vector<double> vL;
    vector<int> binHist;
    for (int i = 0; i <= 35; i++) {
      binHist.push_back(0);
    }
    vector<vector<int>> vHist;
    for (int y = 0; y < im.height(); y++) {
      HSLAPixel* p = im.getPixel(x, y);
      if (x-1 < 0 && y-1 < 0) {
        vHX.push_back(cos(p->h*PI/180));
        vHY.push_back(sin(p->h*PI/180));
        vS.push_back(p->s);
        vL.push_back(p->l);
      } else if (x-1 < 0) {
        vHX.push_back(vHX[y-1]+cos(p->h*PI/180));
        vHY.push_back(vHY[y-1]+sin(p->h*PI/180));
        vS.push_back(vS[y-1]+p->s);
        vL.push_back(vL[y-1]+p->l);
        for (int i = 0; i <= 35; i++) {
          binHist[i] = vHist[y-1][i];
        }
      } else if (y-1 < 0) {
        vHX.push_back(sumHueX[x-1][y]+cos(p->h*PI/180));
        vHY.push_back(sumHueY[x-1][y]+sin(p->h*PI/180));
        vS.push_back(sumSat[x-1][y]+p->s);
        vL.push_back(sumLum[x-1][y]+p->l);
        for (int i = 0; i <= 35; i++) {
          binHist[i] = hist[x-1][y][i];
        }
      } else {
        vHX.push_back(sumHueX[x-1][y]+vHX[y-1]+cos(p->h*PI/180)-sumHueX[x-1][y-1]);
        vHY.push_back(sumHueY[x-1][y]+vHY[y-1]+sin(p->h*PI/180)-sumHueY[x-1][y-1]);
        vS.push_back(sumSat[x-1][y]+vS[y-1]+p->s-sumSat[x-1][y-1]);
        vL.push_back(sumLum[x-1][y]+vL[y-1]+p->l-sumLum[x-1][y-1]);
        binHist[(int)(p->h/10)]++;
        for (int i = 0; i<= 35; i++) {
          binHist[i] = hist[x-1][y][i]+vHist[y-1][i]-hist[x-1][y-1][i];
        }
      }
      binHist[(int)(p->h/10)]++;
      vHist.push_back(binHist);
    }
    sumHueX.push_back(vHX);
    sumHueY.push_back(vHY);
    sumSat.push_back(vS);
    sumLum.push_back(vL);
    hist.push_back(vHist);
  }
}

long stats::rectArea(pair<int,int> ul, pair<int,int> lr) {
  return (lr.first-ul.first+1)*(lr.second-ul.second+1);
}

HSLAPixel stats::getAvg(pair<int,int> ul, pair<int,int> lr) {
  long nPs = rectArea(ul, lr);
  double hX;
  double hY;
  double s;
  double l;
  if (ul.first-1 < 0 && ul.second-1 < 0) {
    hX = sumHueX[lr.first][lr.second]/nPs;
    hY = sumHueY[lr.first][lr.second]/nPs;
    s = sumSat[lr.first][lr.second]/nPs;
    l = sumLum[lr.first][lr.second]/nPs;
  } else if (ul.first-1 < 0) {
    hX = (sumHueX[lr.first][lr.second]-sumHueX[lr.first][ul.second-1])/nPs;
    hY = (sumHueY[lr.first][lr.second]-sumHueY[lr.first][ul.second-1])/nPs;
    s = (sumSat[lr.first][lr.second]-sumSat[lr.first][ul.second-1])/nPs;
    l = (sumLum[lr.first][lr.second]-sumLum[lr.first][ul.second-1])/nPs;
  } else if (ul.second-1 < 0) {
    hX = (sumHueX[lr.first][lr.second]-sumHueX[ul.first-1][lr.second])/nPs;
    hY = (sumHueY[lr.first][lr.second]-sumHueY[ul.first-1][lr.second])/nPs;
    s = (sumSat[lr.first][lr.second]-sumSat[ul.first-1][lr.second])/nPs;
    l = (sumLum[lr.first][lr.second]-sumLum[ul.first-1][lr.second])/nPs;
  } else {
    hX = (sumHueX[lr.first][lr.second]-sumHueX[ul.first-1][lr.second]-
          sumHueX[lr.first][ul.second-1]+sumHueX[ul.first-1][ul.second-1])/nPs;
    hY = (sumHueY[lr.first][lr.second]-sumHueY[ul.first-1][lr.second]-
          sumHueY[lr.first][ul.second-1]+sumHueY[ul.first-1][ul.second-1])/nPs;
    s = (sumSat[lr.first][lr.second]-sumSat[ul.first-1][lr.second]-
          sumSat[lr.first][ul.second-1]+sumSat[ul.first-1][ul.second-1])/nPs;
    l = (sumLum[lr.first][lr.second]-sumLum[ul.first-1][lr.second]-
          sumLum[lr.first][ul.second-1]+sumLum[ul.first-1][ul.second-1])/nPs;
  }
  double h = atan2(hY, hX) * 180 / PI;
  return HSLAPixel(h, s, l);
}

vector<int> stats::buildHist(pair<int,int> ul, pair<int,int> lr) {
  vector<int> selected = hist[lr.first][lr.second];
  if (ul.first-1 < 0 && ul.second-1 < 0) {
    // do nothing
  } else if (ul.first-1 < 0) { // (0, y), nothing to its left
    for (int i = 0; i <= 35; i++) {
      selected[i] = selected[i]-hist[lr.first][ul.second-1][i];
    }
  } else if (ul.second-1 < 0) { // (x, 0), nothing above
    for (int i = 0; i <= 35; i++) {
      selected[i] = selected[i]-hist[ul.first-1][lr.second][i];
    }
  } else {
    for (int i = 0; i <= 35; i++) {
      selected[i] = selected[i]-hist[lr.first][ul.second-1][i]
                    -hist[ul.first-1][lr.second][i]+hist[ul.first-1][ul.second-1][i];
    }
  }
  return selected;
}

// takes a distribution and returns entropy
// partially implemented so as to avoid rounding issues.
double stats::entropy(vector<int> & distn,int area){
    double entropy = 0.;
    for (int i = 0; i < 36; i++) {
        if (distn[i] > 0 )
            entropy += ((double) distn[i]/(double) area)
                                    * log2((double) distn[i]/(double) area);
    }
    return  -1 * entropy;
}

double stats::entropy(pair<int,int> ul, pair<int,int> lr){
  int area = (int) rectArea(ul, lr);
  vector<int> distn = buildHist(ul, lr);
  return entropy(distn, area);
}
