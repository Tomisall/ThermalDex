import tkinter as tk
from tkinter import ttk

class ChemicalInventoryApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Chemical Inventory System")

        # Database
        self.database = []

        # Create notebook
        self.notebook = ttk.Notebook(root)

        # Database Entry Tab
        self.create_database_entry_tab()

        # Molecule View Tab
        self.create_molecule_view_tab()

        self.notebook.pack(expand=1, fill="both")

    def create_database_entry_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text='Database Entry')

        # Entry widgets
        lbl_name = ttk.Label(tab, text='Chemical Name:')
        lbl_formula = ttk.Label(tab, text='Chemical Formula:')
        entry_name = ttk.Entry(tab)
        entry_formula = ttk.Entry(tab)

        # Add to database function
        def add_to_database():
            name = entry_name.get()
            formula = entry_formula.get()

            if name and formula:
                self.database.append({'Name': name, 'Formula': formula})
                entry_name.delete(0, 'end')
                entry_formula.delete(0, 'end')
                self.update_molecule_view()

        # Button to add to database
        btn_add = ttk.Button(tab, text='Add to Database', command=add_to_database)

        # Grid layout
        lbl_name.grid(row=0, column=0, padx=10, pady=10, sticky='w')
        entry_name.grid(row=0, column=1, padx=10, pady=10)
        lbl_formula.grid(row=1, column=0, padx=10, pady=10, sticky='w')
        entry_formula.grid(row=1, column=1, padx=10, pady=10)
        btn_add.grid(row=2, column=0, columnspan=2, pady=10)

    def create_molecule_view_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text='Molecule View')

        # Listbox to display molecules
        listbox_molecules = tk.Listbox(tab)
        listbox_molecules.pack(expand=True, fill='both')

        def update_listbox():
            listbox_molecules.delete(0, 'end')
            for entry in self.database:
                listbox_molecules.insert('end', f"{entry['Name']} ({entry['Formula']})")

        # Update listbox when switching to this tab
        tab.bind("<Visibility>", lambda event: update_listbox())

    def update_molecule_view(self):
        # Update listbox in the Molecule View tab
        self.notebook.tab(1, state='normal')  # Enable the Molecule View tab
        self.notebook.select(1)  # Switch to the Molecule View tab
        self.notebook.tab(1, state='hidden')  # Disable the Molecule View tab again

if __name__ == "__main__":
    root = tk.Tk()
    app = ChemicalInventoryApp(root)
    root.mainloop()
