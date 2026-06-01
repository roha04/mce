#pragma once

#include <array>
#include <cstddef>
#include <vector>

constexpr int kHex20NodeCount = 20;
constexpr int kElementDofCount = 60;
constexpr int kGaussPointCount3D = 27;
constexpr int kGaussPointCount2D = 9;

// Грань η = +1 (верхня Y = Ly), нумерація як у методичці / MeshRenderer.
constexpr int kHexFaceEtaPlus = 5;

struct Node
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Element
{
    int nodes[kHex20NodeCount]{};
};

struct MaterialProperties
{
    double youngModulus = 2.0e11;
    double poissonRatio = 0.3;
};

struct GaussPoint3D
{
    double xi = 0.0;
    double eta = 0.0;
    double zeta = 0.0;
    double weight = 0.0;
};

using ElementStiffnessMatrix = std::array<std::array<double, kElementDofCount>, kElementDofCount>;
using ElementLoadVector = std::array<double, kElementDofCount>;

class Mesh
{
public:
    struct ZPEntry
    {
        int elementIndex = 0;
        int faceIndex = kHexFaceEtaPlus;
        double pressure = 0.0;
    };

    void generateRectangularParallelepiped(double Lx, double Ly, double Lz, int Nx, int Ny, int Nz);

    // ZU (30), ZP (31) — після generateRectangularParallelepiped.
    void buildBoundaryData(double pressure);

    int computeHalfBandwidth() const;

    const std::vector<Node>& getNodes() const { return nodes_; }
    const std::vector<Element>& getElements() const { return elements_; }

    // AKT (7): [0..nqp-1]=X, [nqp..2nqp-1]=Y, [2nqp..3nqp-1]=Z.
    const std::vector<double>& getAKT() const { return akt_; }
    // NT (21): NT[i][e] — глобальний номер i-го локального вузла e-го СЕ.
    const std::vector<std::array<int, kHex20NodeCount>>& getNT() const { return nt_; }
    const std::vector<int>& getZU() const { return zu_; }
    const std::vector<ZPEntry>& getZP() const { return zp_; }

    std::size_t nodeCount() const { return nodes_.size(); }
    std::size_t elementCount() const { return elements_.size(); }

    double getLx() const { return lx_; }
    double getLy() const { return ly_; }
    double getLz() const { return lz_; }
    int getNx() const { return nx_; }
    int getNy() const { return ny_; }
    int getNz() const { return nz_; }

    void clear();

private:
    static int gridNodeIndex(int ix, int iy, int iz, int Nx, int Ny);
    void buildAKTandNT();

    std::vector<Node> nodes_;
    std::vector<Element> elements_;
    std::vector<double> akt_;
    std::vector<std::array<int, kHex20NodeCount>> nt_;
    std::vector<int> zu_;
    std::vector<ZPEntry> zp_;
    double lx_ = 0.0;
    double ly_ = 0.0;
    double lz_ = 0.0;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
};

std::vector<GaussPoint3D> generateGaussPoints3D();
std::vector<GaussPoint3D> generateGaussPoints2D();

void computeHex20ShapeFunctions(double xi, double eta, double zeta, double N[kHex20NodeCount]);
void computeHex20ShapeDerivativesLocal(double xi, double eta, double zeta,
                                       double dN[kHex20NodeCount][3]);

// DFIABG (39): ∂ψ/∂(ξ,η,ζ) у 27 вузлах Гауса.
class GaussShapeDerivativeCache
{
public:
    static const GaussShapeDerivativeCache& instance();

    const std::vector<GaussPoint3D>& gaussPoints() const { return points_; }
    const double* dNLocal(int gp, int node) const;

private:
    GaussShapeDerivativeCache();

    std::vector<GaussPoint3D> points_;
    double dN_[kGaussPointCount3D][kHex20NodeCount][3]{};
};

using DFIABGCache = GaussShapeDerivativeCache;

// DPSITE (46): ψ та ∂ψ/∂h, ∂ψ/∂t на грані η=+1 у 9 вузлах Гауса.
class FaceGaussShapeCache
{
public:
    static const FaceGaussShapeCache& instance();

    const std::vector<GaussPoint3D>& gaussPoints2D() const { return points2d_; }
    double shapeValue(int gp, int node) const;
    double shapeDerivH(int gp, int node) const;
    double shapeDerivT(int gp, int node) const;

private:
    FaceGaussShapeCache();

    std::vector<GaussPoint3D> points2d_;
    double psi_[kGaussPointCount2D][kHex20NodeCount]{};
    double dPsiDh_[kGaussPointCount2D][kHex20NodeCount]{};
    double dPsiDt_[kGaussPointCount2D][kHex20NodeCount]{};
};

using DPSITECache = FaceGaussShapeCache;

struct Jacobian3x3
{
    double data[3][3]{};
    double determinant = 0.0;

    bool invert(double inverse[3][3]) const;
};

Jacobian3x3 computeJacobian(const double dNLocal[kHex20NodeCount][3],
                            const Node elementNodes[kHex20NodeCount]);

void transformShapeDerivativesToGlobal(const double dNLocal[kHex20NodeCount][3],
                                       const double invJ[3][3],
                                       double dNGlobal[kHex20NodeCount][3]);

void buildIsotropicElasticityMatrix6x6(const MaterialProperties& material, double D[6][6]);

ElementStiffnessMatrix assembleElementStiffnessMatrix(
    const Node elementNodes[kHex20NodeCount],
    const MaterialProperties& material,
    const GaussShapeDerivativeCache& gaussCache = GaussShapeDerivativeCache::instance());

ElementLoadVector assembleElementLoadVectorTopFace(const Node elementNodes[kHex20NodeCount],
                                                   double tractionY,
                                                   int faceIndex = kHexFaceEtaPlus);

bool isElementOnTopFace(const Node elementNodes[kHex20NodeCount], double ly, double tolerance = 1.0e-9);
bool isElementOnBottomFace(const Node elementNodes[kHex20NodeCount], double tolerance = 1.0e-9);

void getHex20LocalNodeCoordinates(int localNode, double& xi, double& eta, double& zeta);

// Маркер налагодження (заняття 13): det(J)=8 на одиничному кубі.
double verifyUnitCubeJacobianDeterminant();
