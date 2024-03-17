import random

def generate_matrix_file():
    # Demande à l'utilisateur la taille de la matrice
    n = int(input("Entrez la taille de la matrice : "))

    # Génère une matrice carrée de taille n x n avec des nombres réels entre -1 et 1
    matrix = [[round(random.uniform(-2, 2), 3) for _ in range(n)] for _ in range(n)]

    # Construit le nom du fichier en fonction de la taille de la matrice
    file_name = f"matrice_{n}.txt"

    # Écrit la taille de la matrice au début du fichier
    with open(file_name, "w") as file:
        file.write(str(n) + "\n")

        # Écrit les éléments de la matrice dans le fichier
        for row in matrix:
            file.write(" ".join(map(str, row)) + "\n")

    print(f"Matrice générée et enregistrée dans le fichier : {file_name}")

# Appelle la fonction pour générer la matrice en demandant à l'utilisateur la taille
generate_matrix_file()
