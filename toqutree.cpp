
/**
 *
 * toqutree (pa3)
 * significant modification of a quadtree .
 * toqutree.cpp
 * This file will be used for grading.
 *
 */

#include "toqutree.h"

toqutree::Node::Node(pair<int,int> ctr, int dim, HSLAPixel a)
	:center(ctr),dimension(dim),avg(a),NW(NULL),NE(NULL),SE(NULL),SW(NULL)
	{}

toqutree::~toqutree(){
	clear(root);
	numberOfNodes = 0;
}

toqutree::toqutree(const toqutree & other) {
	root = copy(other.root);
	numberOfNodes = other.numberOfNodes;
}


toqutree & toqutree::operator=(const toqutree & rhs){
	if (this != &rhs) {
		clear(root);
		root = copy(rhs.root);
		numberOfNodes = rhs.numberOfNodes;
	}
	return *this;
}

/*
* This constructor grabs the 2^k x 2^k sub-image centered
* in imIn and uses it to build a quadtree. It may assume
* that imIn is large enough to contain an image of that size.
*/
toqutree::toqutree(PNG & imIn, int k) {
	int ulf = imIn.width()/2-pow(2, k-1);
	int uls = imIn.height()/2-pow(2, k-1);
	int lrf = imIn.width()/2-1+pow(2, k-1);
	int lrs = imIn.height()/2-1+pow(2, k-1);
	if (k-1 < 0) { // i.e. k == 0, 1x1 square
		ulf = 0;
		uls = 0;
		lrf = 0;
		lrs = 0;
	}
	PNG* nPNG = new PNG(pow(2, k), pow(2, k));
	for (int x = ulf; x <= lrf; x++) {
		for (int y = uls; y <= lrs; y++) {
			*nPNG->getPixel(x-ulf, y-uls) = *imIn.getPixel(x, y);
		}
	}
	numberOfNodes = 0;
	root = buildTree(nPNG, k);
}

double toqutree::avgEntropy(PNG* imIn, int k, int x, int y, stats* stats) {
	// calculate entropy of SE square
	int southE = calcEntropy(imIn, k, x, y, stats);
	// calculate entropy of SW square (right)
	int southW = calcEntropy(imIn, k, x+pow(2, k), y, stats);
	// calculate entropy of NE square (below)
	int northE = calcEntropy(imIn, k, x, y+pow(2, k), stats);
	// calculaet entropy of NW square (below right)
	int northW = calcEntropy(imIn, k, x+pow(2, k), y+pow(2, k), stats);
	return (southE+southW+northE+northW)/4;
}

double toqutree::calcEntropy(PNG* im, int k, int x, int y, stats* stats) {
	if (x >= (int)im->width()) {
		x = mod(x, im->width());
	}
	if (y >= (int)im->height()) {
		y = mod(y, im->height());
	}
	int lrf = x-1+pow(2, k);
	int lrs = y-1+pow(2, k);
	auto ul = pair<int, int> (x, y);
	auto lr = pair<int, int> (lrf, lrs);
	if (lrf < (int)im->width() && lrs < (int)im->height()) {
		return stats->entropy(ul, lr);
	} else if (lrf >= (int)im->width() && lrs >= (int)im->height()) { // extra on all four corners
		int extraX = lrf-im->width();
		int extraY = lrs-im->height();
		lr = pair<int, int> ((int)im->width()-1, (int)im->height()-1);
		auto ul1 = pair<int, int> (0, y);
		auto lr1 = pair<int, int> (extraX, (int)im->height()-1);
		auto ul2 = pair<int, int> (0, 0);
		auto lr2 = pair<int, int> (extraX, extraY);
		auto ul3 = pair<int, int> (x, 0);
		auto lr3 = pair<int, int> (im->width()-1, extraY);
		vector<int> total = stats->buildHist(ul, lr);
		vector<int> v2 = stats->buildHist(ul1, lr1);
		vector<int> v3 = stats->buildHist(ul2, lr2);
		vector<int> v4 = stats->buildHist(ul3, lr3);
		for (int i = 0; i <= 35; i++) {
		 total[i] = total[i]+v2[i]+v3[i]+v4[i];
		}
		int area = (int) stats->rectArea(ul, lr)+stats->rectArea(ul1, lr1)+stats->rectArea(ul2, lr2)+stats->rectArea(ul3, lr3);
		return stats->entropy(total, area);
	} else if (lrf >= (int)im->width()) { // extra on left
		int extraX = lrf-im->width();
		lr = pair<int, int> (im->width()-1, lrs);
		auto ul1 = pair<int, int> (0, y);
		auto lr1 = pair<int, int> (extraX, lrs);
		vector<int> total = stats->buildHist(ul, lr);
		vector<int> v2 = stats->buildHist(ul1, lr1);
		for (int i = 0; i <= 35; i++) {
		 total[i] = total[i]+v2[i];
		}
		int area = (int) stats->rectArea(ul, lr) + stats->rectArea(ul1, lr1);
		return stats->entropy(total, area);
	} else {															// extra on top
		int extraY = lrs-im->height();
		lr = pair<int, int> (lrf, (int)im->height()-1);
		auto ul1 = pair<int, int> (x, 0);
		auto lr1 = pair<int, int> (lrf, extraY);
		vector<int> total = stats->buildHist(ul, lr);
		vector<int> v2 = stats->buildHist(ul1, lr1);
		for (int i = 0; i <= 35; i++) {
		 total[i] = total[i]+v2[i];
		}
		int area = (int) stats->rectArea(ul, lr) + stats->rectArea(ul1, lr1);
		return stats->entropy(total, area);
	}
	return 0.0;
}

PNG* toqutree::createSubImage(PNG * im, pair<int, int> ul, int k) {
	PNG* nPNG = new PNG(pow(2, k), pow(2, k));
	// check if ul is out of bounds
	if (ul.first >= (int)im->width()) {
		ul.first = mod(ul.first, im->width());
	}
	if (ul.second >= (int)im->height()) {
		ul.second = mod(ul.second, im->height());
	}

	// calculate lower right corner coordinates
	int lrf = ul.first+pow(2, k)-1;
	int lrs = ul.second+pow(2, k)-1;
	// check if wrapping needed (and wrap) and copy correct pixels over to nPNG
	if (lrf < (int)im->width() && lrs < (int)im->height()) {
		for (int x = ul.first; x <= lrf; x++) {
			for (int y = ul.second; y <= lrs; y++) {
				*nPNG->getPixel(x-ul.first, y-ul.second) = *im->getPixel(x, y);
			}
		}
	} else if (lrf >= (int)im->width() && lrs >= (int)im->height()) { // wrap all four corners
		int leftLRX = lrf-im->width();
		int topLRY = lrs-im->height();
		int w = im->width()-ul.first; // width on nPNG filled by lower right corner
		int h = im->height()-ul.second;	// height on nPNG filled by lower right corner
		auto lr1 = pair<int, int> (im->width()-1, im->height()-1);	// lower right corner
		for (int x = ul.first; x <= lr1.first; x++) {
			for (int y = ul.second; y <= lr1.second; y++) {
				*nPNG->getPixel(x-ul.first, y-ul.second) = *im->getPixel(x, y);
			}
		}
		auto ul2 = pair<int, int> (0, ul.second);										// lower left corner
		auto lr2 = pair<int, int> (leftLRX, im->height()-1);
		for (int x = ul2.first; x <= lr2.first; x++) {
			for (int y = ul2.second; y <= lr2.second; y++) {
				*nPNG->getPixel(x-ul2.first+w, y-ul2.second) = *im->getPixel(x, y);
			}
		}
		auto ul4 = pair<int, int> (ul.first, 0);										// upper right corner
		auto lr4 = pair<int, int> (im->width()-1, topLRY);
		for (int x = ul4.first; x <= lr4.first; x++) {
			for (int y = ul4.second; y <= lr4.second; y++) {
				*nPNG->getPixel(x-ul4.first, y-ul4.second+h) = *im->getPixel(x, y);
			}
		}
		auto ul3 = pair<int, int> (0, 0);										// upper left corner
		auto lr3 = pair<int, int> (leftLRX, topLRY);
		for (int x = ul3.first; x <= lr3.first; x++) {
			for (int y = ul3.second; y <= lr3.second; y++) {
				*nPNG->getPixel(x-ul3.first+w, y-ul3.second+h) = *im->getPixel(x, y);
			}
		}
	} else if (lrf >= (int)im->width()) { // wrap extra on left
		int leftLRX = lrf-im->width();
		lrf = (int)im->width()-1;
		int w = im->width()-ul.first;
		for (int x = ul.first; x <= lrf; x++) {
			for (int y = ul.second; y <= lrs; y++) {
				*nPNG->getPixel(x-ul.first, y-ul.second) = *im->getPixel(x, y);
			}
		}
		auto ul1 = pair<int, int> (0, ul.second);
		auto lr1 = pair<int, int> (leftLRX, lrs);
		for (int x = ul1.first; x <= lr1.first; x++) {
			for (int y = ul1.second; y <= lr1.second; y++) {
				*nPNG->getPixel(x-ul1.first+w, y-ul1.second) = *im->getPixel(x, y);
			}
		}
	} else if (lrs >= (int)im->height()) {	// wrap extra on top
		int topLRY = lrs-im->height();
		lrs = (int)im->height()-1;
		int h = im->height()-ul.second;
		for (int x = ul.first; x <= lrf; x++) {
			for (int y = ul.second; y <= lrs; y++) {
				*nPNG->getPixel(x-ul.first, y-ul.second) = *im->getPixel(x, y);
			}
		}
		auto ul1 = pair<int, int> (ul.first, 0);
		auto lr1 = pair<int, int> (lrf, topLRY);
		for (int x = ul1.first; x <= lr1.first; x++) {
			for (int y = ul1.second; y <= lr1.second; y++) {
				*nPNG->getPixel(x-ul1.first, y-ul1.second+h) = *im->getPixel(x, y);
			}
		}
	}
	return nPNG;
}


// Note that you will want to practice careful memory use
// In this function. We pass the dynamically allocated image
// via pointer so that it may be released after it is used.
// similarly, at each level of the tree you will want to
// declare a dynamically allocated stats object, and free it
// once you've used it to choose a split point, and calculate
// an average.
toqutree::Node * toqutree::buildTree(PNG * im, int k) {
	// set up
	stats* statsObject = new stats(*im);

	// calculate split point
	int ulf = (im->width()/2)-pow(2, k-2);
	int uls = (im->height()/2)-pow(2, k-2);
	int lrf = (im->width()/2)-1+pow(2, k-2);
	int lrs = (im->height()/2)-1+pow(2, k-2);
	if (k <= 0) {
		ulf = 0;
		uls = 0;
		lrf = 0;
		lrs = 0;
	}
	auto sp = pair<int, int>(ulf, uls);	// ul (x, y) of possible split point square
	double lowest = 0.0;
	for (int x = ulf; x <= lrf; x++) {
		for (int y = uls; y <= lrs; y++) {
			// calculate average entropy of four resulting squares
			// if avg entropy at curr is lower than lowest
			// change splitting point to current point
			double curr = avgEntropy(im, k-1, x, y, statsObject);
			if (lowest == 0.0 || curr < lowest) {
				lowest = curr;
				sp = pair<int, int>(x, y);
			}
		}
	}

	// create root of toqutree
	auto ul = pair<int, int>(0, 0);
	auto lr = pair<int, int>(im->width()-1, im->height()-1);
	HSLAPixel a = statsObject->getAvg(ul, lr);		// calculate avg of currenmt entire im
	if (a.h < 0) a.h = a.h+360;
	Node* n = new Node(sp, k, a);
	numberOfNodes++;
	// delete statsObject for proper memory management because we no longer need it

	delete statsObject;
	statsObject = NULL;

	// if current node is 2x2 or bigger, recursively create and assign children
	// if 1x1 (k==0), no need to make + assign children
	if (k >= 1) {
		PNG* sEPNG = createSubImage(im, sp, k-1);
		n->SE = buildTree(sEPNG, k-1);

		auto sW = pair<int, int> (sp.first+pow(2, k-1), sp.second);
		PNG* sWPNG = createSubImage(im, sW, k-1);
		n->SW = buildTree(sWPNG, k-1);

		auto nE = pair<int, int> (sp.first, sp.second+pow(2, k-1));
		PNG* nEPNG = createSubImage(im, nE, k-1);
		n->NE = buildTree(nEPNG, k-1);

		auto nW = pair<int, int> (sp.first+pow(2, k-1), sp.second+pow(2, k-1));
		PNG* nWPNG = createSubImage(im, nW, k-1);
		n->NW = buildTree(nWPNG, k-1);
	}

	// delete im for proper memory management
	delete im;
	im = NULL;

	// return pointer to node made
	return n;
}

int toqutree::mod(int coordinate, int p) {
    int m = coordinate%p;
    return m<0 ? m+p : m;
}

HSLAPixel toqutree::findPixel(Node* node, int x, int y) {
	// Base case: if current node is a leaf, 1x1;
	if (node->dimension == 0) {
		return node->avg;
	}
	// If current node is not a leaf, 2x2 or greater;
	if (node->dimension > 0) {
		int childL = (int)pow(2, node->dimension-1);
		int parentL = (int)pow(2, node->dimension);
		// determine coordinates of children
		auto seUL = node->center;
		auto swUL = pair<int, int>(mod((seUL.first+childL), parentL), seUL.second);
		auto neUL = pair<int, int>(seUL.first, mod((seUL.second+childL), parentL));
		auto nwUL = pair<int, int>(mod((seUL.first+childL), parentL), mod((seUL.second+childL), parentL));
		// find child containing x, y
		// change x, y to corresponding coordinates on child
		// call findPixel recursively on child
		int a;
		int b;
		if (seUL.first < swUL.first || (int)swUL.first == 0) {
			if (seUL.second < neUL.second || (int)neUL.second == 0) {
				// SE child is fully on screen
				if (x >= seUL.first && x < seUL.first+childL) {
					if (y >= seUL.second && y < seUL.second+childL) {
						// pixel in SE
						a = mod(x-seUL.first, parentL);
						b = mod(y-seUL.second, parentL);
						if (node->SE != NULL) return findPixel(node->SE, a, b);
					} else {
						// pixel in NE
						a = mod(x-neUL.first, parentL);
						b = mod(y-neUL.second, parentL);
						if (node->NE != NULL) return findPixel(node->NE, a, b);
					}
				} else {
					if (y >= swUL.second && y < swUL.second+childL) {
						// pixel in SW
						a = mod(x-swUL.first, parentL);
						b = mod(y-swUL.second, parentL);
						if (node->SW != NULL) return findPixel(node->SW, a, b);
					} else {
						// pixel in NW
						a = mod(x-nwUL.first, parentL);
						b = mod(y-nwUL.second, parentL);
						if (node->NW != NULL) return findPixel(node->NW, a, b);
					}
				}
			} else if (seUL.second > neUL.second) {
				// SE child wrapped on top only
				// NE child fully on screen
				if (x >= neUL.first && x < neUL.first+childL) {
					if (y >= neUL.second && y < neUL.second+childL) {
						// pixel in NE
						a = mod(x-neUL.first, parentL);
						b = mod(y-neUL.second, parentL);
						if (node->NE != NULL) return findPixel(node->NE, a, b);
					} else {
						// pixel in SE
						a = mod(x-seUL.first, parentL);
						b = mod(y-seUL.second, parentL);
						if (node->SE != NULL) return findPixel(node->SE, a, b);
					}
				} else {
					if (y >= nwUL.second && y < nwUL.second+childL) {
						// pixel in NW
						a = mod(x-nwUL.first, parentL);
						b = mod(y-nwUL.second, parentL);
						if (node->SE != NULL) return findPixel(node->NW, a, b);
					} else {
						// pixel in SW
						a = mod(x-swUL.first, parentL);
						b = mod(y-swUL.second, parentL);
						if (node->SW != NULL) return findPixel(node->SW, a, b);
					}
				}
			}
		} else if (seUL.first > swUL.first) {
			if (seUL.second < neUL.second || neUL.second == 0) {
				// SE child wrapped on left only
				// SW child fully on screen
				if (x >= swUL.first && x < swUL.first+childL) {
					if (y >= swUL.second && y < swUL.second+childL) {
						// pixel in SW
						a = mod(x-swUL.first, parentL);
						b = mod(y-swUL.second, parentL);
						if (node->SW != NULL) return findPixel(node->SW, a, b);
					} else {
						// pixel in NW
						a = mod(x-nwUL.first, parentL);
						b = mod(y-nwUL.second, parentL);
						if (node->NW != NULL) return findPixel(node->NW, a, b);
					}
				} else {
					if (y >= seUL.second && y < seUL.second+childL) {
						// pixel in SE
						a = mod(x-seUL.first, parentL);
						b = mod(y-seUL.second, parentL);
						if (node->SE != NULL) return findPixel(node->SE, a, b);
					} else {
						// pixel in NE
						a = mod(x-neUL.first, parentL);
						b = mod(y-neUL.second, parentL);
						if (node->NE != NULL) return findPixel(node->NE, a, b);
					}
				}
			} else if (seUL.second > neUL.second){
				// SE child wrapped on all four corners
				// NW child fully on screen
				if (x >= nwUL.first && x < nwUL.first+childL) {
					if (y >= nwUL.second && y < nwUL.second+childL) {
						// pixel in NW
						a = mod(x-nwUL.first, parentL);
						b = mod(y-nwUL.second, parentL);
						if (node->NW != NULL) return findPixel(node->NW, a, b);
					} else {
						// pixel in SW
						a = mod(x-swUL.first, parentL);
						b = mod(y-swUL.second, parentL);
						if (node->SW != NULL) return findPixel(node->SW, a, b);
					}
				} else {
					if (y >= neUL.second && y < neUL.second+childL) {
						// pixel in NE
						a = mod(x-neUL.first, parentL);
						b = mod(y-neUL.second, parentL);
						if (node->NE != NULL) return findPixel(node->NE, a, b);
					} else {
						// pixel in SE
						a = mod(x-seUL.first, parentL);
						b = mod(y-seUL.second, parentL);
						if (node->SE != NULL) return findPixel(node->SE, a, b);
					}
				}
			}
		}
	}
	return node->avg;
}

// My algorithm for this problem included a helper function
// that was analogous to Find in a BST, but it navigated the
// quadtree, instead.
PNG toqutree::render() {
	int side = pow(2, root->dimension);
	PNG nPNG = PNG(side, side);
	for (int x = 0; x < side; x++) {
		for (int y = 0; y < side; y++) {
			*nPNG.getPixel(x, y) = findPixel(root, x, y);
		}
	}
	return nPNG;
}

bool toqutree::withinTol(Node* & node, HSLAPixel & avg, double tol) {
	// if node is greater than 2x2, recursive call
	// if node is 2x2, its children are leaves
	bool se;
	bool sw;
	bool ne;
	bool nw;
	if (node->dimension > 1) {
		se = withinTol(node->SE, avg, tol);
		sw = withinTol(node->SW, avg, tol);
		ne = withinTol(node->NE, avg, tol);
		nw = withinTol(node->NW, avg, tol);
	} else {
		se = node->SE->avg.dist(avg) <= tol;
		sw = node->SW->avg.dist(avg) <= tol;
		ne = node->NE->avg.dist(avg) <= tol;
		nw = node->NW->avg.dist(avg) <= tol;
	}
	return (se && sw && ne && nw);
}


void toqutree::prune(Node* & node, double tol) {
	if (withinTol(node, node->avg, tol)) {
		// all of the nodes children and their leaves are within tolerance
		clear(node->SE);
		clear(node->SW);
		clear(node->NE);
		clear(node->NW);
		cout << "Node pruned." << endl;
	} else if (node->dimension > 1){
		prune(node->SE, tol);
		prune(node->SW, tol);
		prune(node->NE, tol);
		prune(node->NW, tol);
	}
}

void toqutree::prune(double tol) {
	prune(root, tol);
}

int toqutree::size() {
	return numberOfNodes;
}

/* called by destructor and assignment operator*/
void toqutree::clear(Node * & curr) {
	if (curr != NULL) {
		// if node is greater than 2x2 (k == 1), recursively call clear on children
		// else node is 2x2 and its children are leaf nodes. delete leaf nodes.
		// delete current node
		if (curr->dimension > 1) {
			clear(curr->NW);
			clear(curr->NE);
			clear(curr->SE);
			clear(curr->SW);
		} else {
			delete curr->NW;
			delete curr->NE;
			delete curr->SE;
			delete curr->SW;
			curr->NW = NULL;
			curr->NE = NULL;
			curr->SE = NULL;
			curr->SW = NULL;
		}
		delete curr;
		curr = NULL;
	}
}

/* done */
/* called by assignment operator and copy constructor */
toqutree::Node * toqutree::copy(const Node * other) {
	Node* ncopy = new Node(other->center, other->dimension, other->avg);
	if (other->dimension > 1) {
		ncopy->NW = copy(other->NW);
		ncopy->NE = copy(other->NE);
		ncopy->SE = copy(other->SE);
		ncopy->SW = copy(other->SW);
	}
	if (other->dimension == 1) {
		ncopy->NW = new Node(other->NW->center, other->NW->dimension, other->NW->avg);
		ncopy->NE = new Node(other->NE->center, other->NE->dimension, other->NE->avg);
		ncopy->SE = new Node(other->SE->center, other->SE->dimension, other->SE->avg);
		ncopy->SW = new Node(other->SW->center, other->SW->dimension, other->SW->avg);
	}
	return ncopy;
}
