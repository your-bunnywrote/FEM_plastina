#include "Mesh.h"

int main() {
	Mesh mesh;
	mesh.output("nodes.txt");
	mesh.output("elements.txt", mesh.elements);
	return 0;
}