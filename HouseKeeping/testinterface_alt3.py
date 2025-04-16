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

        # Search Database Tab
        self.create_search_database_tab()

        # Add to Database Tab
        self.create_add_to_database_tab()

        self.notebook.pack(expand=1, fill="both")

    def create_search_database_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text='Search Database')

        # Entry widget for searching
        lbl_search = ttk.Label(tab, text='Search:')
        entry_search = ttk.Entry(tab)
        lbl_search.grid(row=0, column=0, padx=10, pady=10, sticky='w')
        entry_search.grid(row=0, column=1, padx=10, pady=10, sticky='we')

        # Listbox to display search results
        listbox_search_results = tk.Listbox(tab)
        listbox_search_results.grid(row=1, column=0, columnspan=2, padx=10, pady=10, sticky='nsew')

        def search_database():
            query = entry_search.get().lower()
            results = [entry for entry in self.database if query in entry['Name'].lower() or query in entry['Formula'].lower()]

            listbox_search_results.delete(0, 'end')
            for result in results:
                listbox_search_results.insert('end', f"{result['Name']} ({result['Formula']})")

        # Button to search
        btn_search = ttk.Button(tab, text='Search', command=search_database)
        btn_search.grid(row=2, column=0, columnspan=2, pady=10)

    def create_add_to_database_tab(self):
        tab = ttk.Frame(self.notebook)
        self.notebook.add(tab, text='Add to Database')

        # Entry widgets for adding to the database
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

        # Button to add to database
        btn_add = ttk.Button(tab, text='Add to Database', command=add_to_database)

        # Grid layout
        lbl_name.grid(row=0, column=0, padx=10, pady=10, sticky='w')
        entry_name.grid(row=0, column=1, padx=10, pady=10)
        lbl_formula.grid(row=1, column=0, padx=10, pady=10, sticky='w')
        entry_formula.grid(row=1, column=1, padx=10, pady=10)
        btn_add.grid(row=2, column=0, columnspan=2, pady=10)

if __name__ == "__main__":
    root = tk.Tk()
    app = ChemicalInventoryApp(root)
    root.mainloop()
