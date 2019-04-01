
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
}

toqutree::toqutree(const toqutree & other) {
	root = copy(other.root);
}


toqutree & toqutree::operator=(const toqutree & rhs){
	if (this != &rhs) {
		clear(root);
		root = copy(rhs.root);
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
	PNG* nPNG = new PNG(pow(2, k), pow(2, k));
	for (int x = ulf; x <= lrf; x++) {
		for (int y = uls; y <= lrs; y++) {
			*nPNG->getPixel(x-ulf, y-uls) = *imIn.getPixel(x, y);
		}
	}
	root = buildTree(nPNG, k);
}

double toqutree::avgEntropy(PNG * imIn, int k, int x, int y, stats * stats) {
	// calculate entropy of SE square
	int southE = calcEntropy(imIn, k, x, y, stats);
	// calculate entropy of SW square (right)
	int southW = calcEntropy(imIn, k, x+pow(2, k)-1, y, stats);
	// calculate entropy of NE square (below)
	int northE = calcEntropy(imIn, k, x, y+pow(2, k)-1, stats);
	// calculaet entropy of NW square (below right)
	int northW = calcEntropy(imIn, k, x+pow(2, k)-1, y+pow(2, k)-1, stats);
	return (southE+southW+northE+northW)/4;
}

double toqutree::calcEntropy(PNG * im, int k, int x, int y, stats * stats) {
	if (x >= (int)im->width()) {
		x = x%im->width();
	}
	if (y >= (int)im->height()) {
		y = y%im->height();
	}
	auto ul = pair<int, int> (x, y);
	auto lr = pair<int, int> (0, 0);
	int lrf = x-1+pow(2, k);
	int lrs = y-1+pow(2, k);
	if (lrf < (int)im->width() && lrs < (int)im->height()) {
		lr = pair<int, int> (lrf, lrs);
		return stats->entropy(ul, lr);
	} else if (lrf >= (int)im->width() && lrs >= (int)im->height()) { // extra on all four corners
		int extraX = lrf-im->width();
		int extraY = lrs-im->height();
		auto ul1 = pair<int, int> (0, y);
		auto lr1 = pair<int, int> (extraX, im->height()-1);
		auto ul2 = pair<int, int> (0, 0);
		auto lr2 = pair<int, int> (extraX, extraY);
		auto ul3 = pair<int, int> (x, 0);
		auto lr3 = pair<int, int> (im->width()-1, extraY);
		lr = pair<int, int> (im->width()-1, im->height()-1);
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
		lr = pair<int, int> (lrf, im->height()-1);
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
}

PNG* toqutree::createSubImage(PNG * im, pair<int, int> ul, int k) {
	PNG* nPNG = new PNG(pow(2, k), pow(2, k));
	// check if ul is out of bounds
	if (ul.first >= (int)im->width()) {
		ul.first = ul.first%im->width();
	}
	if (ul.second >= (int)im->height()) {
		ul.second = ul.second%im->height();
	}
	// cout<< "ul.first:" << ul.first << endl;
	// cout<< "im->width():" << im->width() << endl;
	// cout<< "ul.second:" << ul.second << endl;
	// cout<< "im->height():" << im->height() << endl;

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
		int extraX = lrf-im->width();
		int extraY = lrs-im->height();
		int w = im->width()-1-ul.first; // width on nPNG filled by lower right corner
		int h = im->height()-1-ul.second;	// height on nPNG filled by lower right corner
		auto lr1 = pair<int, int> (im->width()-1, im->height()-1);	// lower right corner
		for (int x = ul.first; x <= lr1.first; x++) {
			for (int y = ul.second; y <= lr1.second; y++) {
				*nPNG->getPixel(x-ul.first, y-ul.second) = *im->getPixel(x, y);
			}
		}
		auto ul2 = pair<int, int> (0, ul.second);										// lower left corner
		auto lr2 = pair<int, int> (extraX, im->height()-1);
		for (int x = ul2.first; x <= lr2.first; x++) {
			for (int y = ul2.second; y <= lr2.second; y++) {
				*nPNG->getPixel(x-ul2.first+w, y-ul2.second) = *im->getPixel(x, y);
			}
		}
		auto ul3 = pair<int, int> (0, 0);										// upper left corner
		auto lr3 = pair<int, int> (extraX, extraY);
		for (int x = ul3.first; x <= lr3.first; x++) {
			for (int y = ul3.second; y <= lr3.second; y++) {
				*nPNG->getPixel(x-ul3.first+w, y-ul3.second+h) = *im->getPixel(x, y);
			}
		}
		auto ul4 = pair<int, int> (ul.first, 0);										// upper right corner
		auto lr4 = pair<int, int> (im->width()-1, extraY);
		for (int x = ul4.first; x <= lr4.first; x++) {
			for (int y = ul4.second; y <= lr4.second; y++) {
				*nPNG->getPixel(x-ul4.first, y-ul4.second+h) = *im->getPixel(x, y);
			}
		}
	} else if (lrf >= (int)im->width()) { // wrap extra on left
		int extraX = lrf-im->width();
		lrf = im->width()-1;
		int w = lrf-ul.first;
		auto ul1 = pair<int, int> (0, ul.second);
		auto lr1 = pair<int, int> (extraX, lrs);
		for (int x = ul.first; x <= lrf; x++) {
			for (int y = ul.second; y <= lrs; y++) {
				*nPNG->getPixel(x-ul.first, y-ul.second) = *im->getPixel(x, y);
			}
		}
		for (int x = ul1.first; x <= lr1.first; x++) {
			for (int y = ul1.second; y <= lr1.second; y++) {
				*nPNG->getPixel(x-ul1.first+w, y-ul1.second) = *im->getPixel(x, y);
			}
		}
	} else {												// wrap extra on top
		int extraY = lrs-im->height();
		lrs = im->height()-1;
		int h = lrs-ul.second;
		auto ul1 = pair<int, int> (ul.first, 0);
		auto lr1 = pair<int, int> (lrf, extraY);
		for (int x = ul.first; x <= lrf; x++) {
			for (int y = ul.second; y <= lrs; y++) {
				*nPNG->getPixel(x-ul.first, y-ul.second) = *im->getPixel(x, y);
			}
		}
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
	auto sp = pair<int, int>(0, 0);
	stats* statsObject = new stats(*im);

	// calculate split point
	int ulf = (im->width()/2)-pow(2, k-2);
	int uls = (im->height()/2)-pow(2, k-2);
	int lrf = (im->width()/2)-1+pow(2, k-2);
	int lrs = (im->height()/2)-1+pow(2, k-2);
	auto ul = pair<int, int>(ulf, uls);	// ul (x, y) of possible split point square
	auto lr = pair<int, int>(lrf, lrs); // lr (x, y) of possible split point square
	double lowest = 0.0;
	for (int x = ulf; x <= lrf; x++) {
		for (int y = uls; y <= lrs; y++) {
			// calculate average entropy of four resulting squares
			// if avg entropy at curr is lower than lowest
			// change splitting point to current point
			double curr = avgEntropy(im, k-1, x, y, statsObject);
			if (lowest == 0.0 || curr < lowest) {
				sp = pair<int, int>(x, y);
			}
		}
	}
	// cout<< "Split point calculated." << endl;

	// create root of toqutree
	ul = pair<int, int>(0, 0);
	lr = pair<int, int>(im->width()-1, im->height()-1);
	HSLAPixel a = statsObject->getAvg(ul, lr);		// calculate avg of current entire im
	Node* n = new Node(sp, k, a);
	numberOfNodes++;
	// cout<< "Node created." << endl;

	// if current node is bigger than 2x2, recursively create and assign children
	if (k > 1) {
		PNG* sEPNG = createSubImage(im, sp, k-1);
		n->SE = buildTree(sEPNG, k-1);
		// cout<< "SE built." << endl;

		auto sW = pair<int, int> (sp.first+pow(2, k-1)-1, sp.second);
		PNG* sWPNG = createSubImage(im, sW, k-1);
		n->SW = buildTree(sWPNG, k-1);
		// cout<< "SW built." << endl;

		auto nE = pair<int, int> (sp.first, sp.second+pow(2, k-1)-1);
		PNG* nEPNG = createSubImage(im, nE, k-1);
		n->NE = buildTree(nEPNG, k-1);
		// cout<< "NE built." << endl;

		auto nW = pair<int, int> (sp.first+pow(2, k-1)-1, sp.second+pow(2, k-1)-1);
		PNG* nWPNG = createSubImage(im, nW, k-1);
		n->NW = buildTree(nWPNG, k-1);
		// cout<< "NW built." << endl;
	}
	// cout<< "Current node is not bigger than 2x2." << endl;

	// if current node is 2x2, its children are leaf nodes
	if (k == 1) {
		auto ctr = pair<int, int>(0, 0); // each node will only have a single pixel at (0, 0)
		// cout << im->width() << endl;
		// cout << im->height() << endl;
		n->SE = new Node(ctr, 0, *im->getPixel(1, 1));
		n->SW = new Node(ctr, 0, *im->getPixel(0, 1));
		n->NE = new Node(ctr, 0, *im->getPixel(1, 0));
		n->NW = new Node(ctr, 0, *im->getPixel(0, 0));
		numberOfNodes = numberOfNodes+4;
	}
	// cout<< "Leaf nodes are assigned." << endl;

	// delete input png for proper memory management
	delete im;

	// return pointer to node made
	return n;
}

PNG* toqutree::render(Node* node) {

	PNG* nPNG = new PNG(pow(2, node->dimension), pow(2, node->dimension));

	// Base case: if current node is a leaf;
	if (node->dimension == 0) {
		*nPNG->getPixel(0, 0) = node->avg;
	}

	// If current node is not a leaf;
	if (node->dimension > 0) {
		PNG* sEPNG = render(node->SE);
		auto ul1 = node->center;
		nPNG = putBack(nPNG, sEPNG, ul1);
		delete sEPNG;

		PNG* sWPNG = render(node->SW);
		auto ul2 = pair<int, int> (ul1.first+pow(2, node->dimension-1)-1, ul1.second);
		nPNG = putBack(nPNG, sWPNG, ul2);
		delete sWPNG;

		PNG* nEPNG = render(node->NE);
		auto ul3 = pair<int, int> (ul1.first, ul1.second+pow(2, node->dimension-1)-1);
		nPNG = putBack(nPNG, nEPNG, ul3);
		delete nEPNG;

		PNG* nWPNG = render(node->NW);
		auto ul4 = pair<int, int> (ul1.first+pow(2, node->dimension-1)-1, ul1.second+pow(2, node->dimension-1)-1);
		nPNG = putBack(nPNG, nWPNG, ul4);
		delete nWPNG;
	}
	return nPNG;
}

PNG* toqutree::putBack(PNG* im, PNG* subImage, pair<int, int> ul) {
	if (ul.first >= (int)im->width()) {
		ul.first = ul.first%im->width();
	}
	if (ul.second >= (int)im->height()) {
		ul.second = ul.second%im->height();
	}
	// calculate lower right corner coordinates
	int lrf = ul.first+subImage->width()-1;
	int lrs = ul.second+subImage->height()-1;
	// check if wrapping needed (and wrap) and copy correct pixels over to im
	if (lrf < (int)im->width() && lrs < (int)im->height()) {
		for (int x = ul.first; x <= lrf; x++) {
			for (int y = ul.second; y <= lrs; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul.first, y-ul.second);
				cout << "nx: " << x << " ny: " << y << endl;
			}
		}
	} else if (lrf >= (int)im->width() && lrs >= (int)im->height()) { // wrap all four corners
		int extraX = lrf-im->width();
		int extraY = lrs-im->height();
		int w = im->width()-1-ul.first;
		int h = im->height()-1-ul.second;
		auto lr1 = pair<int, int> (im->width()-1, im->height()-1);	// lower right corner
		for (int x = ul.first; x <= lr1.first; x++) {
			for (int y = ul.second; y <= lr1.second; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul.first, y-ul.second);
				cout << "bx1: " << x << " by1: " << y << endl;
			}
		}
		auto ul2 = pair<int, int> (0, ul.second);										// lower left corner
		auto lr2 = pair<int, int> (extraX, im->height()-1);
		for (int x = ul2.first; x <= lr2.first; x++) {
			for (int y = ul2.second; y <= lr2.second; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul2.first+w, y-ul2.second);
				cout << "bx2: " << x << " by2: " << y << endl;
			}
		}
		auto ul3 = pair<int, int> (0, 0);										// upper left corner
		auto lr3 = pair<int, int> (extraX, extraY);
		for (int x = ul3.first; x <= lr3.first; x++) {
			for (int y = ul3.second; y <= lr3.second; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul3.first+w, y-ul3.second+h);
				cout << "bx3: " << x << " by3: " << y << endl;
			}
		}
		auto ul4 = pair<int, int> (ul.first, 0);										// upper right corner
		auto lr4 = pair<int, int> (im->width()-1, extraY);
		for (int x = ul4.first; x <= lr4.first; x++) {
			for (int y = ul4.second; y <= lr4.second; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul4.first, y-ul4.second+h);
				cout << "bx4: " << x << " by4: " << y << endl;
			}
		}
	} else if (lrf >= (int)im->width()) { // wrap extra on left
		int extraX = lrf-im->width();
		lrf = im->width()-1;
		int w = lrf-ul.first;
		auto ul1 = pair<int, int> (0, ul.second);
		auto lr1 = pair<int, int> (extraX, lrs);
		for (int x = ul.first; x <= lrf; x++) {
			for (int y = ul.second; y <= lrs; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul.first, y-ul.second);
				cout << "wx1: " << x << " wy1: " << y << endl;
			}
		}
		for (int x = ul1.first; x <= lr1.first; x++) {
			for (int y = ul1.second; y <= lr1.second; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul1.first+w, y-ul1.second);
				cout << "wx2: " << x << " wy2: " << y << endl;
			}
		}
	} else {												// wrap extra on top
		int extraY = lrs-im->height();
		lrs = im->height()-1;
		int h = lrs-ul.second;
		auto ul1 = pair<int, int> (ul.first, 0);
		auto lr1 = pair<int, int> (lrf, extraY);
		for (int x = ul.first; x <= lrf; x++) {
			for (int y = ul.second; y <= lrs; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul.first, y-ul.second);
				cout << "hx1: " << x << " hy1: " << y << endl;
			}
		}
		for (int x = ul1.first; x <= lr1.first; x++) {
			for (int y = ul1.second; y <= lr1.second; y++) {
				*im->getPixel(x, y) = *subImage->getPixel(x-ul1.first, y-ul1.second+h);
				cout << "hx2: " << x << " hy2: " << y << endl;
			}
		}
	}
	return im;
}

// My algorithm for this problem included a helper function
// that was analogous to Find in a BST, but it navigated the
// quadtree, instead.
PNG toqutree::render() {
	return *render(root);
}

/* oops, i left the implementation of this one in the file! */
//	prune(root,tol);
void toqutree::prune(double tol) {

}

int toqutree::size() {
	return numberOfNodes;
}

/* called by destructor and assignment operator*/
void toqutree::clear(Node * & curr) {
/* your code here */
}

/* done */
/* called by assignment operator and copy constructor */
toqutree::Node * toqutree::copy(const Node * other) {
/* your code here */
}
