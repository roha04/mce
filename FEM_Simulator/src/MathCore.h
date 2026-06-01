#pragma once

#include <array>
#include <cstddef>
#include <vector>

constexpr int kHex20NodeCount = 20;
constexpr int kElementDofCount = 60; // 20 вузлів × 3 DOF (u, v, w)
constexpr int kGaussPointCount3D = 27;

// Глобальний вузол сітки: координати в декартовій системі (X, Y, Z).
struct Node
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

// 20-вузловий serendipity-шестигранник (HEX20):
//   індекси 0..7  — кутові вузли;
//   індекси 8..19 — вузли на серединах ребер.
// nodes[i] зберігає глобальний номер i-го локального вузла (A_{i,e}).
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

class Mesh
{
public:
    void generateRectangularParallelepiped(double Lx, double Ly, double Lz,
                                         int Nx, int Ny, int Nz);

    // L = max_e (N_max - N_min + 1) * 3 — півширина стрічки для DOF u,v,w.
    int computeHalfBandwidth() const;

    const std::vector<Node>& getNodes() const { return nodes_; }
    const std::vector<Element>& getElements() const { return elements_; }

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

    std::vector<Node> nodes_;
    std::vector<Element> elements_;
    double lx_ = 0.0;
    double ly_ = 0.0;
    double lz_ = 0.0;
    int nx_ = 0;
    int ny_ = 0;
    int nz_ = 0;
};

// --- Етап 2: локальні координати та інтегрування ---

std::vector<GaussPoint3D> generateGaussPoints3D();

void computeHex20ShapeFunctions(double xi, double eta, double zeta, double N[kHex20NodeCount]);

void computeHex20ShapeDerivativesLocal(double xi, double eta, double zeta,
                                       double dN[kHex20NodeCount][3]);

class GaussShapeDerivativeCache
{
public:
    static const GaussShapeDerivativeCache& instance();

    const std::vector<GaussPoint3D>& gaussPoints() const { return points_; }

    // ∂N_node/∂(ξ,η,ζ) у вузлі Гауса gp.
    const double* dNLocal(int gp, int node) const;

private:
    GaussShapeDerivativeCache();

    std::vector<GaussPoint3D> points_;
    double dN_[kGaussPointCount3D][kHex20NodeCount][3]{};
};

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

using ElementLoadVector = std::array<double, kElementDofCount>;

constexpr int kGaussPointCount2D = 9;

std::vector<GaussPoint3D> generateGaussPoints2D();

// F^e для грані η = +1 (верхня грань Y = Ly): тяга tractionY [Пa] у напрямку +Y.
ElementLoadVector assembleElementLoadVectorTopFace(const Node elementNodes[kHex20NodeCount],
                                                   double tractionY);

bool isElementOnTopFace(const Node elementNodes[kHex20NodeCount], double ly, double tolerance = 1.0e-9);
bool isElementOnBottomFace(const Node elementNodes[kHex20NodeCount], double tolerance = 1.0e-9);

// Локальні координати (ξ,η,ζ) вузла HEX20 у [-1,1]^3.
void getHex20LocalNodeCoordinates(int localNode, double& xi, double& eta, double& zeta);
